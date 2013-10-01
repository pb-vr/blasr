/*
 * ============================================================================
 *
 *       Filename:  HDFCCSReader_gtest.cpp
 *
 *    Description:  Test common/data/hdf/HDFCCSReader.h
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

#include "data/hdf/HDFCCSReader.h"
#include "CCSSequence.h"
#include "gtest/gtest.h"
#include "../../testdata.h"

using namespace std;
using namespace H5;

class HDFCCSReaderTEST : public ::testing::Test {
public:
    virtual void SetUp() {
    }
    virtual void TearDown() {
    }
};

TEST_F(HDFCCSReaderTEST, ReadCCSFromBasH5) {
    string fileName = baxFile2;
    HDFCCSReader<CCSSequence> reader;
    reader.InitializeDefaultIncludedFields();
    ASSERT_EQ(reader.Initialize(fileName), 1);
    ASSERT_EQ(reader.GetMovieName(),
            "m130220_114643_42129_c100471902550000001823071906131347_s1_p0");

    CCSSequence seq;
    for(int i=0; i < 1000; i++) {
        reader.GetNext(seq);
    }
    reader.Close();
}

TEST_F(HDFCCSReaderTEST, ReadCCSFromCCSH5) {
    string fileName = ccsFile1;
    HDFCCSReader<CCSSequence> reader;
    reader.SetReadBasesFromCCS();
    reader.InitializeDefaultIncludedFields();
    ASSERT_EQ(reader.Initialize(fileName), 1);
    ASSERT_EQ(reader.GetMovieName(),
            "m130328_211423_ethan_c100499512550000001823070408081371_s1_p0");

    CCSSequence seq;
    for(int i=0; i < 1000; i++) {
        reader.GetNext(seq);
    }

    reader.Close();
}



