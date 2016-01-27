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

#include "IConverter.h"
#include <pbbam/BamRecord.h>
#include <algorithm>
#include <iostream>
#include <set>
#include <cassert>
#include <cmath>
using namespace PacBio;
using namespace PacBio::BAM;
using namespace std;

IConverter::IConverter(Settings& settings)
    : settings_(settings)
{ }

IConverter::~IConverter(void) { }

void IConverter::AddErrorMessage(const std::string& e)
{ errors_.push_back(e); }

BamHeader IConverter::CreateHeader(const string& modeString)
{
    BamHeader header;

    // @HD VN:<current SAM/BAM spec version>
    //     SO:unsorted
    //     pb:<current PacBio BAM spec version>
    header.Version("1.5")
          .SortOrder("unknown")
          .PacBioBamVersion("3.0.1");

    // @RG ID: <read group ID>
    //     DS: READTYPE=<HQREGION|POLYMERASE|SUBREAD>[;<Tag Manifest>;BINDINGKIT=<foo>;SEQUENCINGKIT=<bar>;BASECALLERVERSION=<42>]
    //     PL: PACBIO
    //     PU: <movieName>
    //
    ReadGroupInfo rg(settings_.movieName, modeString);
    rg.BindingKit(bindingKit_)
      .SequencingKit(sequencingKit_)
      .BasecallerVersion(basecallerVersion_)
      .FrameRateHz(frameRateHz_);

    if (settings_.usingDeletionQV)      rg.BaseFeatureTag(BaseFeature::DELETION_QV,      "dq");
    if (settings_.usingDeletionTag)     rg.BaseFeatureTag(BaseFeature::DELETION_TAG,     "dt");
    if (settings_.usingInsertionQV)     rg.BaseFeatureTag(BaseFeature::INSERTION_QV,     "iq");
    if (settings_.usingMergeQV)         rg.BaseFeatureTag(BaseFeature::MERGE_QV,         "mq");
    if (settings_.usingSubstitutionQV)  rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV,  "sq");
    if (settings_.usingSubstitutionTag) rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_TAG, "st");
    if (settings_.usingIPD) {
        FrameCodec codec = FrameCodec::V1;
        if (settings_.losslessFrames)
            codec = FrameCodec::RAW;
        rg.IpdCodec(codec, "ip");
    }
    if (settings_.usingPulseWidth) {
        FrameCodec codec = FrameCodec::V1;
        if (settings_.losslessFrames)
            codec = FrameCodec::RAW;
        rg.PulseWidthCodec(codec, "pw");
    }

    header.AddReadGroup(rg);

    // @PG ID:bax2bam-<version>
    //     PN:bax2bam
    //     CL:bax2bam <args>
    //     DS:<description>
    //     VN:<version>

    ProgramInfo program(settings_.program + "-" + settings_.version);
    program.Name(settings_.program)
           .CommandLine(settings_.program + " " + settings_.args)
           .Description(settings_.description)
           .Version(settings_.version);
    header.AddProgram(program);

    return header;
}

std::vector<std::string> IConverter::Errors(void) const
{ return errors_; }
