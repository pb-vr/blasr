/*
 * =====================================================================================
 *
 *       Filename:  SeqUtils_gtest.cpp
 *
 *    Description:  Test common/SeqUtils.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:21:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include "gtest/gtest.h"
#include "SeqUtils.h"

TEST(SeqUtils, OnlyACTG){
    DNASequence seq;
    Nucleotide seqnt[] = "ATGC";
    seq.seq = seqnt;
    seq.length = 4;
    EXPECT_EQ(OnlyACTG(seq), 1);

    Nucleotide seqnt1[] = "ATXYZ";
    seq.seq = seqnt1;
    seq.length = 5;
    EXPECT_EQ(OnlyACTG(seq), 0);
}

TEST(SeqUtils, CountMasked){
    DNASequence seq;
    Nucleotide seqnt[] = "ATGCNNNNNNATGC";
    seq.seq = seqnt;
    seq.length = 14;
    EXPECT_EQ(CountMasked(seq), 6);
}

TEST(SeqUtils, CountNotMasked){
    DNASequence seq;
    Nucleotide seqnt[] = "ATGCNNNNNNATGC";
    seq.seq = seqnt;
    seq.length = 14;
    EXPECT_EQ(CountNotMasked(seq), 8);

}
