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

// Author: Yuan Li

#include "TestData.h"

#include "SMRTSequence.hpp"
#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFZMWWriter.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

TEST(HDFZMWWriter, EndToEnd)
{
    // setup
    std::string outfn = tests::Out_Dir  + "/" + "zmw.h5";

    HDFFile outfile;
    outfile.Open(outfn, H5F_ACC_TRUNC);

    HDFZMWWriter writer(outfn, outfile.rootGroup);

    for (int i = 1 ; i < 1000; i++) {
        SMRTSequence seq;
        seq.length = i;
        seq.zmwData.holeNumber = i;
        seq.zmwData.holeStatus = static_cast<unsigned char> (8);
        seq.readScore = 0.87;
        seq.HQRegionSnr('A', 0.1);
        seq.HQRegionSnr('C', 0.2);
        seq.HQRegionSnr('G', 0.3);
        seq.HQRegionSnr('T', 0.4);

        bool OK = writer.WriteOneZmw(seq);
        EXPECT_TRUE(OK);
    }
};
