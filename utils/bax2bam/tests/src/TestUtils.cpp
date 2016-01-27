// Copyright (c) 2014, Pacific Biosciences of California, Inc.
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

#include "TestData.h"
#include "TestUtils.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

void RemoveFiles(const vector<string>& filenames)
{
    for (auto fn : filenames)
        remove(fn.c_str());
}

void RemoveFile(const string& filename)
{
    vector<string> filenames;
    filenames.push_back(filename);
    RemoveFiles(filenames);
}

int RunBax2Bam(const vector<string>& baxFilenames,
               const string& outputType,
               const string& additionalArgs)
{
    string convertArgs;
    convertArgs += outputType;
    if (!additionalArgs.empty()) {
        convertArgs += string(" ");
        convertArgs += additionalArgs;
    }
    for (auto fn : baxFilenames) {
        convertArgs += string(" ");
        convertArgs += fn;
    }

    const string& convertCommandLine = tests::Bax2Bam_Exe + string(" ") + convertArgs;
    return system(convertCommandLine.c_str());
}
