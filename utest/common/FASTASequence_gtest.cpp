/*
 * =====================================================================================
 *
 *       Filename:  FASTASequence_gtest.cpp
 *
 *    Description:  Test common/FASTASequence.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:19:13 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "FASTASequence.h"
#include "gtest/gtest.h"
#include <climits>
#include <iostream>
#include <fstream>

class FASTASequenceTest : public ::testing::Test {
public:
    virtual void SetUp() {
    }

    virtual void TearDown() {
    }

    FASTASequence fastaOne;
    FASTASequence fastaTwo;
    FASTASequence fastaThree;

    std::streambuf * sbuf;
    ofstream ofs;
};


//Test FASTASequence 
TEST_F(FASTASequenceTest, ALLFUNC) {
    // Test constructor 
    EXPECT_TRUE(fastaOne.title == NULL);
    EXPECT_TRUE(fastaOne.titleLength == 0);
    EXPECT_TRUE(fastaOne.seq == NULL);
    EXPECT_TRUE(fastaOne.length == 0);
    EXPECT_FALSE(fastaOne.deleteOnExit);
    EXPECT_EQ(fastaOne.GetStorageSize(), 0);

    DNASequence dna;
    Nucleotide thisNuc[] = "ATGCATGCTC";
    dna.seq = thisNuc;
    dna.length = 10;

    int titleLength = 22;
    string title("fasta_seq_one comments");
    fastaOne.title = new char [titleLength];
    memcpy(fastaOne.title, title.c_str(), titleLength);

    fastaOne.titleLength = titleLength;

    EXPECT_EQ(fastaOne.GetName(), string("fasta_seq_one"));
    fastaOne.seq = thisNuc;
    fastaOne.length = 10;

    // use ShallowCopy carefully, since title may double free
    // fastaTwo.ShallowCopy(fastaOne);
    // EXPECT_EQ(fastaTwo.seq, fastaOne.seq);
    // EXPECT_EQ(fastaTwo.title, fastaOne.title);


    // Test AppendToTitle
    fastaOne.AppendToTitle(string("XXX"));
    EXPECT_EQ(fastaOne.titleLength, 26);

    string newTitle = "fasta_seq_one commentsXXX";
    EXPECT_STREQ(fastaOne.title, newTitle.c_str());


    // Test ReverseComplementSelf()
    fastaOne.ReverseComplementSelf();
    string rcSeq = "GAGCATGCAT";
    for (int i = 0; i < rcSeq.size(); i++) {
        EXPECT_EQ(fastaOne.seq[i], rcSeq[i]);
    }

    // Test operator =
    fastaTwo=fastaOne;
    EXPECT_NE(fastaOne.title, fastaTwo.title);
    EXPECT_EQ(fastaOne.titleLength, fastaTwo.titleLength);
    EXPECT_STREQ(fastaOne.title, fastaTwo.title);
    EXPECT_NE(fastaOne.seq, fastaTwo.seq);
    for (int i = 0; i < fastaOne.length; i++) {
        EXPECT_EQ(fastaOne.seq[i], fastaTwo.seq[i]);
    }

    // Test MakeRC(rhs&)
    fastaOne.MakeRC(fastaThree);

    // Test PrintSeq
    stringstream ss;
    fastaThree.PrintSeq(ss);
    string thisTitle, thisComment, thisSeq;
    ss >> thisTitle;
    ss >> thisComment;
    ss >> thisSeq;

    EXPECT_EQ(thisTitle, ">fasta_seq_one");
    EXPECT_EQ(thisComment, "commentsXXX");
    EXPECT_EQ(thisSeq, "ATGCATGCTC");
}
