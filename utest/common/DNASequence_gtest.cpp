/*
 * =====================================================================================
 *
 *       Filename:  DNASequence_gtest.cpp
 *
 *    Description:  Test cpp/common/DNASequence.h
 *
 *        Version:  1.0
 *        Created:  10/27/2012 09:42:13 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "DNASequence.h"
#include <climits>
#include <iostream>
#include <fstream>
//using namespace std;

//Note ::testing::Test not ::testing::TEST
//SetUp() and TearDown(), not Setup() and Teardown()
class DNASequenceTest : public ::testing::Test {
public:
    DNASequence dnaOne;
};

//Test DNASequence constructor
TEST_F(DNASequenceTest, Constructor) {
    DNASequence dnaSeq;
    EXPECT_TRUE(dnaSeq.seq == NULL);
    EXPECT_TRUE(dnaSeq.length == 0);
    EXPECT_TRUE(dnaSeq.size() == dnaSeq.length);
    EXPECT_TRUE(dnaSeq.bitsPerNuc == 8);
    EXPECT_FALSE(dnaSeq.deleteOnExit);

    Nucleotide HKITTY[] = "HELLO,KITTY!";
    dnaSeq.seq = HKITTY;
    dnaSeq.length = sizeof(HKITTY)/sizeof(Nucleotide) - 1;
//    dnaSeq.Print(cout);
    EXPECT_EQ(dnaSeq.size(), 12);

    
    DNALength thisLen = 12;
    Nucleotide * thisNuc = new Nucleotide [thisLen];
    memcpy(thisNuc, HKITTY, thisLen);
    DNASequence newDnaSeq; 
    newDnaSeq.seq = thisNuc;
    newDnaSeq.length = thisLen;
//    newDnaSeq.Print(cout);
    EXPECT_EQ(memcmp(newDnaSeq.seq, dnaSeq.seq, thisLen), 0);
    EXPECT_EQ(newDnaSeq.length, thisLen);
    if (!thisNuc) delete thisNuc;

    DNASequence nnewDnaSeq;
    thisLen = 12;
    string atgc ("atgcatgcatgc");
    thisNuc = new Nucleotide [thisLen];
    for (int i = 0 ; i < thisLen; i++) {
        thisNuc[i] = atgc[i];
    }
    string ret;
    nnewDnaSeq.seq = thisNuc;
    nnewDnaSeq.length = thisLen;
    for (int i = 0 ; i < thisLen; i++) {
        ret += nnewDnaSeq.seq[i];
    }
    EXPECT_STREQ(ret.c_str(), atgc.c_str());
}

//Test DNASequence Append()
TEST_F(DNASequenceTest, Append) {
    DNALength oneLen = 10;
    Nucleotide * one = new Nucleotide [oneLen];

    string As("AAAAAAAAAA"); 
    for (int i = 0; i < oneLen; i++) {
        one[i] = As[i];
    }
    //Can not memcpy a string to a DNASequence directly 
    //such as memcpy(one, As.c_str(), oneLen), because 
    //DNASequence.seq is of type unsigned char, not char

    DNALength twoLen = 20;
    Nucleotide * two = new Nucleotide [twoLen];
    
    string Gs("GGGGGGGGGGGGGGGGGGGG");
    for (int i = 0; i < twoLen; i++) {
        two[i] = Gs[i];
    }
    //memcpy(two, Gs.c_str(), twoLen);

    Nucleotide * three = new Nucleotide [oneLen + twoLen];
    memcpy(three, one, oneLen);
    memcpy(three+oneLen, two, twoLen);

    DNASequence dnaTwo;
    dnaOne.seq = one; 
    dnaOne.length = oneLen;
    dnaTwo.seq = two;
    dnaTwo.length = twoLen;

    dnaOne.Append(dnaTwo, 0);
    EXPECT_EQ(dnaOne.length, oneLen + twoLen);
    EXPECT_EQ(memcmp(dnaOne.seq, three, dnaOne.length), 0);

    string AGs("AAAAAAAAAAGGGGGGGGGGGGGGGGGGGG");
    for (int i = 0; i < dnaOne.length; i++) {
        EXPECT_EQ(AGs[i], (char)dnaOne.seq[i]);
    }
  
    //if appendPos is positive, overwrite this sequence
    //from appendPos to the end
    AGs = "AAGGGGGGGGGGGGGGGGGGGG";
    DNALength appendPos = 2;
    dnaOne.Append(dnaTwo, appendPos);
    EXPECT_EQ(dnaOne.length, appendPos + twoLen);
    for (int i = 0; i < dnaOne.length; i++) {
        EXPECT_EQ(AGs[i], (char)dnaOne.seq[i]);
    }
    
    if(!one) delete one;
    if(!two) delete two;
    if(!three) delete three;
}

//Test DNASequence TakeOwnership
TEST_F(DNASequenceTest, TakeOwnership) {
    DNALength oneLen = 10;
    Nucleotide * one = new Nucleotide [oneLen];
    
    dnaOne.seq = one; 
    dnaOne.length = oneLen;

    DNASequence dnaTwo;
    dnaTwo.TakeOwnership(dnaOne);
    EXPECT_EQ(dnaTwo.length, dnaOne.length);
    EXPECT_EQ(dnaTwo.deleteOnExit, dnaOne.deleteOnExit);
    EXPECT_EQ(dnaTwo.seq, dnaOne.seq);

    dnaTwo.deleteOnExit = true;
    dnaTwo.TakeOwnership(dnaOne);
    //a bug may occur if deleteOneExit is true and 
    //TakeOwnership() is called twice. In that case, both
    //dnaOne and dnaTwo will become wild pointers 
    EXPECT_EQ(dnaOne.seq, dnaTwo.seq);
    if(!one) delete one;
}

//Test DNASequence ShallowCopy
TEST_F(DNASequenceTest, ShallowCopy) {
    DNALength oneLen = 10;
    Nucleotide * one = new Nucleotide [oneLen];

    string As("AAAAAAAAAA");
    for (int i = 0; i < oneLen; i++) {
        one[i] = As[i];
    }
    dnaOne.seq = one;
    dnaOne.length = oneLen;

    DNASequence dnaTwo;
    dnaTwo.ShallowCopy(dnaOne);

    EXPECT_EQ(dnaTwo.length, dnaOne.length);
    EXPECT_EQ(dnaTwo.seq   , dnaOne.seq);
    EXPECT_EQ(dnaTwo.deleteOnExit, dnaOne.deleteOnExit);
}


//Test DNASequence.Copy(const DNASequence rhs, 
//                      DNALength rhsPos,
//                      DNALength rhsLength)
TEST_F(DNASequenceTest, Copy) {
    DNALength oneLen = 10;
    Nucleotide * one = new Nucleotide [oneLen];

    string As("AGAAAAACAA");
    for (int i = 0; i < oneLen; i++) {
        one[i] = As[i];
    }

    dnaOne.seq = one;
    dnaOne.length = oneLen;

    DNASequence dnaTwo;
    dnaTwo.Copy(dnaOne);

    EXPECT_EQ(dnaTwo.length, dnaOne.length);
    EXPECT_NE(dnaTwo.seq   , dnaOne.seq);
    EXPECT_TRUE(dnaTwo.deleteOnExit); 
    EXPECT_EQ(memcmp(dnaTwo.seq, dnaOne.seq, dnaOne.length), 0);

    //if rhs.length is 0, return this * 
    DNASequence dnaThree;
    dnaTwo.Copy(dnaThree);
    //dnaTwo remains unchanged
    EXPECT_EQ(dnaTwo.length, 0);
    EXPECT_NE(dnaTwo.seq, dnaOne.seq);
    EXPECT_TRUE(dnaTwo.deleteOnExit); 
    EXPECT_TRUE(dnaTwo.seq == NULL);

    //if rhsPos is not 0 and rhsLength is 0
    dnaTwo.Copy(dnaOne, 2);
    EXPECT_EQ(dnaTwo.length, dnaOne.length - 2);
    EXPECT_TRUE(dnaTwo.deleteOnExit); 
    EXPECT_EQ(memcmp(dnaTwo.seq, dnaOne.seq + 2, dnaTwo.length), 0);


    //if the subsequence to copy is out of bounds
    EXPECT_GT(200, dnaOne.length);
    //EXPECT_EXIT(dnaTwo.Copy(dnaOne, 200), ::testing::ExitedWithCode(1), ""); 


    //if both rhsPos and rhsLength are less than MAXINT,
    //but rhsPos+ rhsLength > MAXINT
    DNALength rhsPos = 3;
    DNALength rhsLength = UINT_MAX -1;
    EXPECT_TRUE(rhsPos < UINT_MAX && rhsLength < UINT_MAX);
    EXPECT_TRUE(rhsLength > dnaOne.length + 1);
    //EXPECT_EXIT(dnaTwo.Copy(dnaOne, rhsPos, rhsLength), ::testing::ExitedWithCode(1), "");


    //if rhsPos > rhs.length
    //EXPECT_EXIT(dnaTwo.Copy(dnaOne, dnaOne.length + 1), ::testing::ExitedWithCode(1), "")
    //    << "Copy a subsequence which is out of bounds. This needs to be taken care of. See bug 21867.";

}

//Test DNASequence Allocate(DNALength)
TEST_F(DNASequenceTest, Allocate) {
    dnaOne.Allocate(0);
    EXPECT_EQ(dnaOne.length, 0);

    DNASequence dnaTwo;
    dnaTwo.Allocate(100);
    EXPECT_EQ(dnaTwo.length, 100);
}

//Test DNASequence ReferenceSubstring(rhs, pos, substrLength)
TEST_F(DNASequenceTest, ReferenceSubstring) {
    DNALength oneLen = 10;
    dnaOne.seq = new Nucleotide[oneLen];
    dnaOne.length = oneLen;

    DNASequence dnaTwo;
    dnaTwo.ReferenceSubstring(dnaOne);

    EXPECT_EQ(dnaOne.seq, dnaTwo.seq);
    EXPECT_EQ(dnaOne.length, dnaTwo.length);
    EXPECT_FALSE(dnaTwo.deleteOnExit);

//    EXPECT_DEATH_IF_SUPPORTED(dnaTwo.ReferenceSubstring(dnaOne, 100), "");
    delete dnaOne.seq;
}

TEST_F(DNASequenceTest, TheNext) {
    EXPECT_TRUE(true);
}

