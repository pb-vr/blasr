// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

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
