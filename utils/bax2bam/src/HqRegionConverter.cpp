// Author: Derek Barnett

#include "HqRegionConverter.h"

#include "utils/RegionUtils.hpp"
#include "HDFRegionTableReader.hpp"
#include <boost/scoped_ptr.hpp>
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <set>
#include <sstream>

using namespace std;
using namespace PacBio::BAM;

HqRegionConverter::HqRegionConverter(Settings& settings)
    : ConverterBase(settings)
{ }

HqRegionConverter::~HqRegionConverter(void) { }

bool HqRegionConverter::ConvertFile(HDFBasReader* reader,
                                    PacBio::BAM::BamWriter* writer)
{
    return ConvertFile(reader, writer, nullptr);
}

bool HqRegionConverter::ConvertFile(HDFBasReader* reader,
                                    PacBio::BAM::BamWriter* writer,
                                    PacBio::BAM::BamWriter* scrapsWriter) 
{
    assert(reader);

    // read region table info
    boost::scoped_ptr<HDFRegionTableReader> regionTableReader(new HDFRegionTableReader);
    RegionTable regionTable;
    std::string fn = filenameForReader_[reader];
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
    int hqStart, hqEnd, score;
    while (reader->GetNext(smrtRecord)) {

        // attempt get high quality region
        if (!LookupHQRegion(smrtRecord.zmwData.holeNumber,
                            regionTable,
                            hqStart,
                            hqEnd,
                            score))
        {
            stringstream s;
            s << "could not find HQ region for hole number: " << smrtRecord.zmwData.holeNumber;
            AddErrorMessage(s.str());
            smrtRecord.Free();
            return false;
        }

        // Catch and repair 1-off errors in the HQ region
        hqEnd = (hqEnd == static_cast<int>(smrtRecord.length)-1) ? smrtRecord.length
                                                                 : hqEnd;

        // sequencing ZMW
        if (IsSequencingZmw(smrtRecord))
        {
            // write HQRegion to main BAM file
            if (hqStart < hqEnd)
            {
                if (!WriteRecord(smrtRecord,
                                 hqStart,
                                 hqEnd,
                                 ReadGroupId(),
                                 writer))
                {
                    smrtRecord.Free();
                    return false;
                }
            }

            // if scraps BAM file present
            if (scrapsWriter)
            {
                // write 5'-end LQ sequence
                if (hqStart > 0)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               0,
                                               hqStart,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write 3'-end LQ sequence
                if (static_cast<size_t>(hqEnd) < smrtRecord.length)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               hqEnd,
                                               smrtRecord.length,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }
            }
        }

        // non-sequencing ZMW
        else
        {
            assert(!IsSequencingZmw(smrtRecord));

            // only write these if scraps BAM present & we are in 'internal mode'
            if (settings_.isInternal && scrapsWriter)
            {
                // write 5'-end LQ sequence
                if (hqStart > 0)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               0,
                                               hqStart,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write HQRegion to scraps BAM file
                if (hqStart < hqEnd)
                {
                    if (!WriteFilteredRecord(smrtRecord,
                                             hqStart,
                                             hqEnd,
                                             ScrapsReadGroupId(),
                                             scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }

                // write 3'-end LQ sequence
                if (static_cast<size_t>(hqEnd) < smrtRecord.length)
                {
                    if (!WriteLowQualityRecord(smrtRecord,
                                               hqEnd,
                                               smrtRecord.length,
                                               ScrapsReadGroupId(),
                                               scrapsWriter))
                    {
                        smrtRecord.Free();
                        return false;
                    }
                }
            }
        }

        smrtRecord.Free();
    }

    // if we get here, all OK
    return true;
}

string HqRegionConverter::HeaderReadType(void) const
{ return "HQREGION"; }

string HqRegionConverter::ScrapsReadType(void) const
{ return "SCRAP"; }

string HqRegionConverter::OutputFileSuffix(void) const
{ return ".hqregions.bam"; }

string HqRegionConverter::ScrapsFileSuffix(void) const
{ return ".lqregions.bam"; }
