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

#include "PolymeraseReadConverter.h"
#include <pbbam/BamRecord.h>
#include <pbbam/BamWriter.h>

using namespace std;

PolymeraseReadConverter::PolymeraseReadConverter(Settings& settings)
    : ConverterBase(settings)
{ }

PolymeraseReadConverter::~PolymeraseReadConverter(void) { }

bool PolymeraseReadConverter::ConvertFile(HDFBasReader* reader,
                                          PacBio::BAM::BamWriter* writer)
{
    assert(reader);

    // initialize BamRecord with default values (shared across all reads)
    PacBio::BAM::BamRecordImpl bamRecord;

    // initialize read scores
    InitReadScores(reader);

    // fetch records from HDF5 file
    SMRTSequence smrtRecord;
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

bool PolymeraseReadConverter::ConvertFile(HDFBasReader* reader,
                                          PacBio::BAM::BamWriter* writer,
                                          PacBio::BAM::BamWriter* scrapsWriter) 
{ return false; }

string PolymeraseReadConverter::HeaderReadType(void) const
{ return "POLYMERASE"; }

string PolymeraseReadConverter::ScrapsReadType(void) const
{ return "UNKNOWN"; }

string PolymeraseReadConverter::OutputFileSuffix(void) const
{ return ".polymerase.bam"; }

string PolymeraseReadConverter::ScrapsFileSuffix(void) const
{ return ".empty.bam"; }
