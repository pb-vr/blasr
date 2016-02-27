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
