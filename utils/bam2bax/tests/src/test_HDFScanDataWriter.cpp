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
#include "TestConstants.h"

#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFScanDataWriter.hpp"
#include "HDFScanDataReader.hpp"
#include "reads/ScanData.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;

const std::string outfn = tests::Out_Dir  + "/" + "scandata.h5";

TEST(HDFScanDataWriter, EndToEnd)
{
    std::map<char, size_t> baseMap = tests::baseMap;

    ScanData scandata;
    tests::make_scandata(scandata, baseMap);
    
    // Write
    tests::write_to(outfn, scandata);

    // Read
    ScanData scandata_;
    EXPECT_TRUE(tests::read_from(outfn, scandata_));

    EXPECT_EQ(scandata_.MovieName(), tests::movieName);
    EXPECT_EQ(scandata_.RunCode(), tests::runCode);
    EXPECT_EQ(scandata_.WhenStarted(), tests::whenStarted);
    EXPECT_EQ(scandata_.FrameRate(), tests::frameRate);
    EXPECT_EQ(scandata_.NumFrames(), tests::numFrames);

    EXPECT_EQ(scandata_.BaseMap()['A'], baseMap['A']);
    /*EXPECT_EQ(scandata_.BaseMap()['C'], tests::baseMap['C']);
    EXPECT_EQ(scandata_.BaseMap()['G'], tests::baseMap['G']);
    EXPECT_EQ(scandata_.BaseMap()['T'], tests::baseMap['T']);*/

    EXPECT_EQ(scandata_.BindingKit(), tests::bindingKit);
    EXPECT_EQ(scandata_.SequencingKit(), tests::sequencingKit);
};

TEST(HDFScanDataWriter, EndToEnd_RandomBaseMap)
{
    std::map<char, size_t> baseMap = tests::randomBaseMap;

    ScanData scandata;
    tests::make_scandata(scandata, baseMap);
    
    // Write
    tests::write_to(outfn, scandata);

    // Read
    ScanData scandata_;
    EXPECT_TRUE(tests::read_from(outfn, scandata_));

    EXPECT_EQ(scandata_.MovieName(), tests::movieName);
    EXPECT_EQ(scandata_.RunCode(), tests::runCode);
    EXPECT_EQ(scandata_.WhenStarted(), tests::whenStarted);
    EXPECT_EQ(scandata_.FrameRate(), tests::frameRate);
    EXPECT_EQ(scandata_.NumFrames(), tests::numFrames);

    EXPECT_EQ(scandata_.BaseMap()['A'], baseMap['A']);
    /*
    EXPECT_EQ(scandata_.BaseMap()['C'], tests::randomBaseMap['C']);
    EXPECT_EQ(scandata_.BaseMap()['G'], tests::randomBaseMap['G']);
    EXPECT_EQ(scandata_.BaseMap()['T'], tests::randomBaseMap['T']);
    */

    EXPECT_EQ(scandata_.BindingKit(), tests::bindingKit);
    EXPECT_EQ(scandata_.SequencingKit(), tests::sequencingKit);
};
