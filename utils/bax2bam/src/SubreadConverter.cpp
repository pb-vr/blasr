// Author: Derek Barnett


#include "SubreadConverter.h"
#include "utils/RegionUtils.hpp"
#include "HDFRegionTableReader.hpp"

#include <boost/scoped_ptr.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

#include <algorithm>

#define MAX( A, B )     ( (A)>(B) ? (A) : (B) )
#define MAX3( A, B, C ) MAX( MAX( A, B ), C )

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

SubreadConverter::SubreadConverter(Settings& settings)
    : ConverterBase(settings)
{ }

SubreadConverter::~SubreadConverter(void) { }

struct SubreadInterval
{
    size_t Start;
    size_t End;
    PacBio::BAM::LocalContextFlags LocalContextFlags;

    SubreadInterval()
        : Start{0}
        , End{0}
        , LocalContextFlags{NO_LOCAL_CONTEXT}
    { }

    SubreadInterval(size_t start, size_t end, bool adapterBefore = false, bool adapterAfter = false)
        : Start{start}
        , End{end}
        , LocalContextFlags{(adapterBefore ? ADAPTER_BEFORE : NO_LOCAL_CONTEXT) |
                            (adapterAfter  ? ADAPTER_AFTER  : NO_LOCAL_CONTEXT)}
    { }
};

inline
bool RegionComparer(const RegionAnnotation& lhs, const RegionAnnotation& rhs)
{
    constexpr int HoleNumber  = RegionAnnotation::HOLENUMBERCOL;
    constexpr int RegionType  = RegionAnnotation::REGIONTYPEINDEXCOL;
    constexpr int RegionStart = RegionAnnotation::REGIONSTARTCOL;

    if (lhs.row[HoleNumber] < rhs.row[HoleNumber])
        return true;
    else if (lhs.row[HoleNumber] == rhs.row[HoleNumber])
    {
        if (lhs.row[RegionType] < rhs.row[RegionType])
            return true;
        else if (lhs.row[RegionType] == rhs.row[RegionType])
            return lhs.row[RegionStart] < rhs.row[RegionStart];
    }
    return false;
}

SubreadInterval ComputeSubreadIntervals(vector<SubreadInterval>* const intervals,
                                        vector<SubreadInterval>* const adapters,
                                        RegionTable& regionTable,
                                        const unsigned holeNumber,
                                        const size_t readLength)
{
    constexpr int RegionStart = RegionAnnotation::REGIONSTARTCOL;
    constexpr int RegionEnd   = RegionAnnotation::REGIONENDCOL;

    // clear the input first
    intervals->clear();

    // region annotations of a zmw
    RegionAnnotations zmwRegions = regionTable[holeNumber];

    // Has non-empty HQregion or not?
    if (not zmwRegions.HasHQRegion())
        return SubreadInterval(0, 0);

    size_t hqStart = zmwRegions.HQStart();
    size_t hqEnd   = zmwRegions.HQEnd();

    // Catch and repair 1-off errors in the HQ region
    hqEnd = (hqEnd == readLength-1) ? readLength : hqEnd;

    // Catch empty or invalid HQ regions and return empty
    if (hqEnd <= hqStart)
        return SubreadInterval(0, 0);

    // adapter intervals of this zmw
    vector<ReadInterval> adapterIntervals = zmwRegions.AdapterIntervals();

    size_t subreadStart  = hqStart;
    bool   adapterBefore = false;

    for (size_t i = 0; i < adapterIntervals.size(); i++) {

        size_t adapterStart = adapterIntervals[i].start;
        size_t adapterEnd   = adapterIntervals[i].end;

        // if we're not in the HQRegion yet, skip ahead
        if (hqStart > adapterEnd)
            continue;

        // if the adapter is beyond the HQRegion, we're done
        if (hqEnd < adapterStart)
            break;

        // If the subread is greater than length=0, save it
        if (subreadStart < adapterStart)
            intervals->emplace_back(SubreadInterval(subreadStart, adapterStart, adapterBefore, true));

        // Save the region of the adapter that overlaps the HQ region
        adapters->emplace_back(SubreadInterval(MAX3(adapterStart, hqStart, subreadStart), 
                    min(adapterEnd, hqEnd)));

        subreadStart  = adapterEnd;
        adapterBefore = true;
    }

    // Save any region between the last adatper and the end of the HQ region as a subread
    if (subreadStart < hqEnd)
        intervals->emplace_back(SubreadInterval(subreadStart, hqEnd, adapterBefore, false));

    return SubreadInterval(hqStart, hqEnd);
}

bool SubreadConverter::ConvertFile(HDFBasReader* reader,
                                   PacBio::BAM::BamWriter* writer)
{
    return ConvertFile(reader, writer, nullptr);
}

bool SubreadConverter::ConvertFile(HDFBasReader* reader,
                                   PacBio::BAM::BamWriter* writer,
                                   PacBio::BAM::BamWriter* scrapsWriter) 
{
    assert(reader);

    // initialize with default values (shared across all unmapped subreads)
    BamRecordImpl bamRecord;

    // read region table info
    boost::scoped_ptr<HDFRegionTableReader> regionTableReader(new HDFRegionTableReader);
    RegionTable regionTable;
    string fn = filenameForReader_[reader];
    assert(!fn.empty());
    if (regionTableReader->Initialize(fn) == 0) {
        AddErrorMessage("could not read region table on "+fn);
        return false;
    }
    regionTable.Reset();
    regionTableReader->ReadTable(regionTable);
    regionTableReader->Close();

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    SMRTSequence smrtRecord;
    while (reader->GetNext(smrtRecord)) {

        SubreadInterval hqInterval;
        vector<SubreadInterval> subreadIntervals;
        vector<SubreadInterval> adapterIntervals;

        // loop over subreads
        try {
            hqInterval = ComputeSubreadIntervals(&subreadIntervals,
                                                 &adapterIntervals,
                                                 regionTable,
                                                 smrtRecord.zmwData.holeNumber,
                                                 smrtRecord.length);
        } catch (runtime_error& e) {
            AddErrorMessage(string(e.what()));
            smrtRecord.Free();
            return false;
        }

        // Write records for the subreads to the primary output BAM
        const size_t endIndex = subreadIntervals.size();
        for (size_t i = 0; i < endIndex; ++i) 
        {
            const int     subreadStart = subreadIntervals[i].Start;
            const int     subreadEnd   = subreadIntervals[i].End;
            const uint8_t contextFlags = subreadIntervals[i].LocalContextFlags;

            // skip invalid or 0-sized intervals
            if (subreadEnd <= subreadStart)
                continue;

            if (!WriteSubreadRecord(smrtRecord, subreadStart, subreadEnd, 
                                    ReadGroupId(), contextFlags, writer))
            {
                smrtRecord.Free();
                return false;
            }
        }

        // Write a record for any 5'-end LQ sequence
        if (scrapsWriter && hqInterval.Start > 0) {
            if (!WriteLowQualityRecord(smrtRecord, 0, hqInterval.Start, ScrapsReadGroupId(), scrapsWriter))
            {
                smrtRecord.Free();
                return false;
            }
        } 

        // Write records for the adapters to the scraps output BAM
        if (scrapsWriter)
        {
            for (size_t i = 0; i < adapterIntervals.size(); ++i) {

                const int adapterStart = adapterIntervals[i].Start;
                const int adapterEnd   = adapterIntervals[i].End;

                // skip invalid or 0-sized adapters
                if (adapterEnd <= adapterStart)
                    continue;

                if (!WriteAdapterRecord(smrtRecord, adapterStart, adapterEnd, ScrapsReadGroupId(), scrapsWriter))
                {
                    smrtRecord.Free();
                    return false;
                }
            }
        }

        // Write a record for any 3'-end LQ sequence
        if (scrapsWriter && hqInterval.End < smrtRecord.length) {
            if (!WriteLowQualityRecord(smrtRecord, hqInterval.End, smrtRecord.length, ScrapsReadGroupId(), scrapsWriter))
            {
                smrtRecord.Free();
                return false;
            }
        } 

        smrtRecord.Free();
    }

    // if we get here, all OK
    return true; 
} 

string SubreadConverter::HeaderReadType(void) const
{ return "SUBREAD"; }

string SubreadConverter::ScrapsReadType(void) const
{ return "SCRAP"; }

string SubreadConverter::OutputFileSuffix(void) const
{ return ".subreads.bam"; }

string SubreadConverter::ScrapsFileSuffix(void) const
{ return ".scraps.bam"; }
