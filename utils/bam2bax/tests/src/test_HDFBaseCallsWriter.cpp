// Author: Yuan Li

#include "TestData.h"
#include "TestConstants.h"

#include "SMRTSequence.hpp"
#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFBaseCallsWriter.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;

const std::string outfn = tests::Out_Dir  + "/" + "basecalls.h5";

TEST(HDFBaseCallsWriter, WriteOneZmw_EndToEnd)
{
    unsigned int len = tests::len;
    // setup a sequence 
    SMRTSequence seq;
    tests::make_smrtseq(seq);
    EXPECT_EQ(seq.length, len);

    // write the sequence to outfn
    bool ret = tests::write_to(outfn, seq);
    EXPECT_TRUE(ret);

    // read the seq from outfn
    SMRTSequence seq2;
    int count = tests::read_from(outfn, seq2);
    EXPECT_EQ(count, 1);

    // compare
    EXPECT_EQ(seq2.length, len);
    EXPECT_EQ(memcmp(seq.seq, seq2.seq, len * sizeof(char)), 0);
    EXPECT_EQ(seq2.zmwData.holeNumber, tests::holeNumber);
    EXPECT_EQ(seq2.zmwData.holeStatus, tests::holeStatus);

    // HQRegionSNR can only be correctly read when ScanData is also available.
    EXPECT_EQ(seq2.HQRegionSnr('A'), -1); 
    EXPECT_EQ(seq2.HQRegionSnr('C'), -1);
    EXPECT_EQ(seq2.HQRegionSnr('G'), -1);
    EXPECT_EQ(seq2.HQRegionSnr('T'), -1);

    EXPECT_EQ(seq2.readScore * 1000, tests::readScoreInt);

    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.deletionQV,  seq2.deletionQV,    len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.deletionTag, seq2.deletionTag,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.insertionQV, seq2.insertionQV,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.mergeQV,     seq2.mergeQV,       len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.substitutionQV,  seq2.substitutionQV,  len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.substitutionTag, seq2.substitutionTag, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.preBaseFrames,    seq2.preBaseFrames, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.widthInFrames,    seq2.widthInFrames, len));
};
