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

    // initialize with default values (shared across all unmapped subreads)
    PacBio::BAM::BamRecordImpl bamRecord;

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
        hqEnd = (hqEnd == static_cast<int>(smrtRecord.length)-1) ? smrtRecord.length : hqEnd;

        // attempt convert BAX to BAM
        if (hqStart < hqEnd)
        {
            // attempt convert BAX to BAM
            if (!WriteRecord(smrtRecord, hqStart, hqEnd, ReadGroupId(), writer))
            {
                smrtRecord.Free();
                return false;
            }
        }

        // Write a record for any 5'-end LQ sequence
        if (scrapsWriter &&  hqStart > 0) {
            if (!WriteLowQualityRecord(smrtRecord, 0, hqStart, ScrapsReadGroupId(), scrapsWriter))
            {
                smrtRecord.Free();
                return false;
            }
        } 

        // Write a record for any 3'-end LQ sequence
        if (scrapsWriter && static_cast<size_t>(hqEnd) < smrtRecord.length) {
            if (!WriteLowQualityRecord(smrtRecord, hqEnd, smrtRecord.length, ScrapsReadGroupId(), scrapsWriter))
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

string HqRegionConverter::HeaderReadType(void) const
{ return "HQREGION"; }

string HqRegionConverter::ScrapsReadType(void) const
{ return "SCRAP"; }

string HqRegionConverter::OutputFileSuffix(void) const
{ return ".hqregions.bam"; }

string HqRegionConverter::ScrapsFileSuffix(void) const
{ return ".lqregions.bam"; }
