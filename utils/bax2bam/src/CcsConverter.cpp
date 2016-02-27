// Author: Derek Barnett

#include <iostream>

#include "CcsConverter.h"

#include "utils/RegionUtils.hpp"
#include "HDFRegionTableReader.hpp"
#include <boost/scoped_ptr.hpp>

#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>
#include <algorithm>

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

CcsConverter::CcsConverter(Settings& settings)
    : ConverterBase(settings)
{
    settings_.usingMergeQV         = false;
    settings_.usingDeletionTag     = false;
    settings_.usingSubstitutionTag = false;
    settings_.usingIPD             = false;
    settings_.usingPulseWidth      = false;
}

CcsConverter::~CcsConverter(void) { }

bool CcsConverter::ConvertFile(HdfCcsReader* reader,
                               BamWriter* writer)
{
    assert(reader);

    // initialize with default values (shared across all unmapped subreads)
    BamRecordImpl bamRecord;

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    CCSSequence smrtRecord;
    while (reader->GetNext(smrtRecord)) {

        // Skip empty records
        if (smrtRecord.length == 0)
            continue;

        // attempt convert BAX to BAM
        if (!WriteRecord(smrtRecord, 0, smrtRecord.length, ReadGroupId(), writer))
        {
            smrtRecord.Free();
            return false;
        }

        smrtRecord.Free();
    }

    // if we get here, all OK
    return true;
}

bool CcsConverter::ConvertFile(HdfCcsReader* reader,
                               PacBio::BAM::BamWriter* writer,
                               PacBio::BAM::BamWriter* scrapsWriter) 
{ return false; }

void CcsConverter::SetSequenceAndQualities(PacBio::BAM::BamRecordImpl* bamRecord,
                                           const CCSSequence& smrtRead,
                                           const int start,
                                           const int length)
{
    recordSequence_.assign((const char*)smrtRead.seq + start, length);
    if (smrtRead.qual.Empty())
        bamRecord->SetSequenceAndQualities(recordSequence_);
    else
    {
        recordQVs_.assign((uint8_t*)smrtRead.qual.data + start,
                          (uint8_t*)smrtRead.qual.data + start + length);
        bamRecord->SetSequenceAndQualities(recordSequence_, recordQVs_.Fastq());
    }
}

void CcsConverter::AddRecordName(PacBio::BAM::BamRecordImpl* bamRecord,
                                 const UInt holeNumber,
                                 const int start,
                                 const int end)
{
    const string name = settings_.movieName + "/"
                      + to_string(holeNumber) + "/ccs";
    bamRecord->Name(name);
}

void CcsConverter::AddModeTags(PacBio::BAM::TagCollection* tags,
                               const CCSSequence& smrtRead,
                               const int start,
                               const int end)
{
    (*tags)["np"] = static_cast<int32_t>(smrtRead.numPasses);
}

CcsConverter::HdfCcsReader* CcsConverter::InitHdfReader()
{
    HdfCcsReader* reader = ConverterBase<CCSSequence, HdfCcsReader>::InitHdfReader();
    // set the reader to CCS mode
    reader->SetReadBasesFromCCS();
    return reader;
}

string CcsConverter::HeaderReadType(void) const
{ return "CCS"; }

string CcsConverter::ScrapsReadType(void) const
{ return "UNKNOWN"; }

string CcsConverter::OutputFileSuffix(void) const
{ return ".ccs.bam"; }

string CcsConverter::ScrapsFileSuffix(void) const
{ return ".empty.bam"; }
