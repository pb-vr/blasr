/*
 * =====================================================================================
 *
 *       Filename:  FASTQReader_gtest.cpp
 *
 *    Description:  Test common/FASTQReader.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:19:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */



#include "gtest/gtest.h"
#include "FASTQReader.h"
#include "testdata.h"

const string movie = 
    "m130328_211423_ethan_c100499512550000001823070408081371_s1_p0";
const int numSeqs = 208;

class FASTQReaderTest:public::testing::Test{
public:
    void SetUp() {
        string filename(fastqFile1);
        reader.Initialize(filename);
    }

    void TearDown() {
        reader.Close();
        seq.Free();
    }

    FASTQReader reader;
    FASTQSequence seq;
};

TEST_F(FASTQReaderTest, GetNext) {
    reader.GetNext(seq);
    EXPECT_EQ(strcmp(seq.title, string(movie + "/8").c_str()), 0);
    EXPECT_EQ(seq.length, 752);
    string expected_seq = string(
        "AATAAAAAAAAAAGAAAGCTTCGAAGTGAGCGAATTACTCTCAGGCAACT"
        "GCGGGTGAAGCCAGAGCAGGCATGATGACACTGGGGAATTTACGCAAATT"
        "TTACCATTGAATTTACACATGCGATGTGCTGGAATGCGGAAGACGGAAAC"
        "GAAACCAGCAATACATCAAACGCCGCACCAGAGAAGAGATATTTGCGCCC"
        "TAAACTAGGTAAGGCGGTTGACTTGAACAGCAAATCAAACGTCAACGAGC"
        "AGCGTGAGTATATACAAGTTATCTCGGATGGAGAACGTATTCTAAATGTA"
        "AGCACGAATCCCGGAAGAGGAAACCAGTTTCTTGGTTTTTCGCCATCCTC"
        "GAAGACCTGTTACAAACCGCACTGGACCTGGAAAGTTTCTGCGCGTAATC"
        "GACAAGACTAGTAACTATCGACATCAACCATCGATTACGGGTTGGGTCAA"
        "TGGGTTCAGATGCAGGTGAGTATCCTTCATATGATAGTCTGACGCTGGCA"
        "TTCGCTCAAAGGAAGTAGACGGTTTTGTAAATAGAAACGCTTGTGAAAAG"
        "CTGAATTTCGCGTCGTCTTCCAGCGATGCAGAGCTGTAGTAGTTCAGATG"
        "ATGACCGTTACTCAAAGTGCCTGCAACGGCTCGGGGCGTGCGCGTCCTGT"
        "GGTGGCTGCTTTTGTTGCGCTGTTTGCAGTGTATGGTTGTCGGGTGATGT"
        "TGCCTGCAAACCCACAAAACCCCACACACACAACAGTTGGGTTGTTGATT"
        "GG");
    string expected_qual = string(
        "(,)'(')''++),.$\"+*$'--.-/+&.$-./$',-.&#'/,.,)-,--,");

    for(int i=0; i < seq.length; i++){
        EXPECT_EQ(seq.seq[i], expected_seq[i]);
    }

    for(int i=0; i < expected_qual.size(); i++){
        EXPECT_EQ(seq.qual[i] + FASTQSequence::charToQuality,
                  expected_qual[i]);
    }

    reader.GetNext(seq);
    EXPECT_EQ(strcmp(seq.title, string(movie+"/9").c_str()), 0);

    // Continue to read
    for(int i=2; i < numSeqs; i++) {
        EXPECT_TRUE(reader.GetNext(seq));
    }
    EXPECT_EQ(strcmp(seq.title, string(movie+"/249").c_str()), 0);

    // Can not proceed.
    EXPECT_FALSE(reader.GetNext(seq));
}

