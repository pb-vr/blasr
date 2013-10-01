/*
 * ============================================================================
 *
 *       Filename:  HDFBasReader_gtest.cpp
 *
 *    Description:  Test common/data/hdf/HDFBasReader.h
 *
 *        Version:  1.0
 *        Created:  08/23/2013 10:17:14 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */

#include "gtest/gtest.h"
#include "data/hdf/HDFBasReader.h"
#include "../../testdata.h"

using namespace std;
using namespace H5;

class HDFBasReaderTEST : public ::testing::Test {
public:
    virtual void SetUp() {
        fileName = baxFile2;
        reader.InitializeDefaultIncludedFields();
        ASSERT_EQ(reader.Initialize(fileName), 1);
    }
    virtual void TearDown() {
        reader.Close();
    }
    string fileName;
    T_HDFBasReader<SMRTSequence> reader; 
};

TEST_F(HDFBasReaderTEST, ReadBaseFromBaseCalls) {
    ASSERT_EQ(reader.GetMovieName(), 
            "m130220_114643_42129_c100471902550000001823071906131347_s1_p0");
    SMRTSequence seq;

    for(int i=0; i < 1000; i++) {
        reader.GetNext(seq);
    }
}


