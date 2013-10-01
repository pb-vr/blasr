/*
 * ============================================================================
 *
 *       Filename:  HDFPlsReader_gtest.cpp
 *
 *    Description:  Test common/data/hdf/HDFPlsReader.h
 *
 *        Version:  1.0
 *        Created:  08/23/2013 11:13:34 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */

#include "data/hdf/HDFPlsReader.h"
#include "gtest/gtest.h"
#include "../../testdata.h"

using namespace std;
using namespace H5;

class HDFPlsReaderTEST : public ::testing::Test {
public:
    virtual void SetUp() {
        fileName = plsFile1;
        ASSERT_EQ(reader.Initialize(fileName), 1);
    }
    virtual void TearDown() {
        reader.Close();
    }
    string fileName;
    HDFPlsReader reader; 
};

TEST_F(HDFPlsReaderTEST, ReadToPulseFile) {
    PulseFile  pulseFile;
    reader.IncludeField("NumEvent");
    reader.IncludeField("StartFrame");
    reader.ReadPulseFileInit(pulseFile);
    reader.ReadPulseFile(pulseFile);
    ASSERT_EQ(pulseFile.platformId, 0);
    ASSERT_EQ(pulseFile.startFrame.size(), 750921871);
}


