/*
 * =====================================================================================
 *
 *       Filename:  FASTAReader_gtest.cpp
 *
 *    Description:  Test common/FASTAReader.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:18:52 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "FASTAReader.h"
#include "testdata.h"

class FASTAReaderTest:public::testing::Test{
public:
    void SetUp() {
        string filename(fastaFile1);
        reader.Initialize(filename);
    }

    void TearDown() {
        reader.Close();
        seq.Free();
    }

    FASTAReader reader;
    FASTASequence seq;
};

TEST_F(FASTAReaderTest, GetNext) {
    reader.GetNext(seq);
    EXPECT_EQ(strcmp(seq.title, "read1"), 0);
    EXPECT_EQ(seq.length, 100);
    string expected_seq = string(
        "AAAAAGGGGGCCCCCACGGCAGCCAGATTTAAATTGAGGGCCCCCCCTTT"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    EXPECT_EQ(strcmp((char*)seq.seq, expected_seq.c_str()), 0); 

    reader.GetNext(seq);
    EXPECT_EQ(strcmp(seq.title, "read2"), 0);
    EXPECT_EQ(seq.length, 99);
}

TEST_F(FASTAReaderTest, ReadAllSequences) {
    vector<FASTASequence> seqs;
    reader.ReadAllSequences(seqs);

    EXPECT_EQ(seqs.size(), 12);
    EXPECT_EQ(strcmp(seqs[11].title, "read2x"), 0);

    string expected_seq = string(
        "AAAAAGGGGGCCCACGGCAGCCAGATTTAAATTGAGGGCAACCCCCCTTT"
        "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
    EXPECT_EQ(strcmp((char*)seqs[11].seq, expected_seq.c_str()), 0);
}

