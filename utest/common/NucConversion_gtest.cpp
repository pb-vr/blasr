/*
 * =====================================================================================
 *
 *       Filename:  NucConversion_gtest.cpp
 *
 *    Description:  Test common/NucConversion.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:21:12 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include "gtest/gtest.h"
#include "NucConversion.h"


// ACGT = 0123
TEST(NucConversion, ASCIITo2BIT) {
    EXPECT_EQ(TwoBit['A'], 0);
    EXPECT_EQ(TwoBit['a'], 0);
    EXPECT_EQ(TwoBit['C'], 1);
    EXPECT_EQ(TwoBit['c'], 1);
    EXPECT_EQ(TwoBit['G'], 2);
    EXPECT_EQ(TwoBit['g'], 2);
    EXPECT_EQ(TwoBit['T'], 3);
    EXPECT_EQ(TwoBit['t'], 3);
    EXPECT_EQ(TwoBit['N'], 255);
    EXPECT_EQ(TwoBit['x'], 255);
}

TEST(NucConversion, ASCIITo3BIT) {
    EXPECT_EQ(ThreeBit['A'], 0);
    EXPECT_EQ(ThreeBit['a'], 0);
    EXPECT_EQ(ThreeBit['C'], 1);
    EXPECT_EQ(ThreeBit['c'], 1);
    EXPECT_EQ(ThreeBit['G'], 2);
    EXPECT_EQ(ThreeBit['g'], 2);
    EXPECT_EQ(ThreeBit['T'], 3);
    EXPECT_EQ(ThreeBit['t'], 3);
    
    EXPECT_EQ(ThreeBit['U'], 4);
    EXPECT_EQ(ThreeBit['M'], 4);
    EXPECT_EQ(ThreeBit['R'], 4);
    EXPECT_EQ(ThreeBit['W'], 4);
    EXPECT_EQ(ThreeBit['S'], 4);
    EXPECT_EQ(ThreeBit['Y'], 4);
    EXPECT_EQ(ThreeBit['K'], 4);
    EXPECT_EQ(ThreeBit['V'], 4);
    EXPECT_EQ(ThreeBit['H'], 4);
    EXPECT_EQ(ThreeBit['D'], 4);
    EXPECT_EQ(ThreeBit['N'], 4);

    EXPECT_EQ(ThreeBit['$'], 5);

    EXPECT_EQ(ThreeBit['p'], 255);
    EXPECT_EQ(ThreeBit['q'], 255);
}

TEST(NucConversion, ISACTG) {
    EXPECT_TRUE(IsACTG['A']);
    EXPECT_TRUE(IsACTG['a']);
    EXPECT_TRUE(IsACTG['C']);
    EXPECT_TRUE(IsACTG['c']);
    EXPECT_TRUE(IsACTG['G']);
    EXPECT_TRUE(IsACTG['g']);
    EXPECT_TRUE(IsACTG['T']);
    EXPECT_TRUE(IsACTG['t']);

    EXPECT_FALSE(IsACTG['w']);
    EXPECT_FALSE(IsACTG['N']);
}


TEST(NucConversion, ThreeBitToASCII) {
    char alphabeta[] = {'A', 'C', 'G', 'T'};
    for(int i = 0; i < 4; i++) {
        EXPECT_EQ(alphabeta[i], ThreeBitToAscii[ThreeBit[alphabeta[i]]]);
    }
}

TEST(NucConversion, AllToUpper) {
    char alphabeta[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
    for(int i = 0; i < 8; i++) {
        EXPECT_EQ(toupper(alphabeta[i]), AllToUpper[alphabeta[i]]);
    }
}

TEST(NucConversion, AllToLower) {
    char alphabeta[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
    for(int i = 0; i < 8; i++) {
        EXPECT_EQ(tolower(alphabeta[i]), AllToLower[alphabeta[i]]);
    }
}

TEST(NucConversion, ReverseComplementNuc) {
    EXPECT_EQ(ReverseComplementNuc['A'], 'T');
    EXPECT_EQ(ReverseComplementNuc['T'], 'A');
    EXPECT_EQ(ReverseComplementNuc['G'], 'C');
    EXPECT_EQ(ReverseComplementNuc['C'], 'G');
}
