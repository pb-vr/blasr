// Author: Yuan Li

#include "TestData.h"
#include "TestConstants.h"

#include "SMRTSequence.hpp"
#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFBasReader.hpp"
#include "HDFBaxWriter.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;

const std::string outfn = tests::Out_Dir  + "/" + "test_HDFBaxWriter.bax.h5";

TEST(HDFBaxWriter, WriteOneZmw_EndToEnd)
{
    // setup a sequence 
    SMRTSequence seq;
    tests::make_smrtseq(seq);
    EXPECT_EQ(seq.length, tests::len);

    // setup a scandata
    ScanData scandata;
    tests::make_scandata(scandata, tests::baseMap);

    // write the sequence to outfn
    bool ret = tests::write_to(outfn, scandata, seq);
    EXPECT_TRUE(ret);

    // read the seq from outfn
    SMRTSequence seq2;
    int count = tests::read_from(outfn, seq2);
    EXPECT_EQ(count, 1);

    // compare
    unsigned int len = tests::len;
    EXPECT_EQ(seq2.length, len);
    EXPECT_EQ(memcmp(seq.seq, seq2.seq, len * sizeof(char)), 0);
    EXPECT_EQ(seq2.zmwData.holeNumber, tests::holeNumber);
    EXPECT_EQ(seq2.zmwData.holeStatus, tests::holeStatus);

    EXPECT_EQ(tests::snra, seq2.HQRegionSnr('A'));
    EXPECT_EQ(tests::snrc, seq2.HQRegionSnr('C'));
    EXPECT_EQ(tests::snrg, seq2.HQRegionSnr('G'));
    EXPECT_EQ(tests::snrt, seq2.HQRegionSnr('T'));

    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.deletionQV,  seq2.deletionQV,    len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.deletionTag, seq2.deletionTag,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.insertionQV, seq2.insertionQV,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.mergeQV,     seq2.mergeQV,       len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.substitutionQV,  seq2.substitutionQV,  len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.substitutionTag, seq2.substitutionTag, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.preBaseFrames,    seq2.preBaseFrames, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.widthInFrames,    seq2.widthInFrames, len));
};

TEST(HDFBaxWriter, WriteOneZmw_EndToEnd_RandomBaseMap)
{
    // setup a sequence 
    SMRTSequence seq;
    tests::make_smrtseq(seq);
    EXPECT_EQ(seq.length, tests::len);

    // setup a scandata
    ScanData scandata;
    tests::make_scandata(scandata, tests::randomBaseMap);

    // write the sequence to outfn
    bool ret = tests::write_to(outfn, scandata, seq);
    EXPECT_TRUE(ret);

    // read the seq from outfn
    SMRTSequence seq2;
    int count = tests::read_from(outfn, seq2);
    EXPECT_EQ(count, 1);

    // compare
    unsigned int len = tests::len;
    EXPECT_EQ(seq2.length, len);
    EXPECT_EQ(memcmp(seq.seq, seq2.seq, len * sizeof(char)), 0);
    EXPECT_EQ(seq2.zmwData.holeNumber, tests::holeNumber);
    EXPECT_EQ(seq2.zmwData.holeStatus, tests::holeStatus);

    EXPECT_EQ(tests::snra, seq2.HQRegionSnr('A'));
    EXPECT_EQ(tests::snrc, seq2.HQRegionSnr('C'));
    EXPECT_EQ(tests::snrg, seq2.HQRegionSnr('G'));
    EXPECT_EQ(tests::snrt, seq2.HQRegionSnr('T'));

    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.deletionQV,  seq2.deletionQV,    len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.deletionTag, seq2.deletionTag,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.insertionQV, seq2.insertionQV,   len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.mergeQV,     seq2.mergeQV,       len));
    EXPECT_TRUE(tests::CmpData<QualityValueVector<unsigned char>>(seq.substitutionQV,  seq2.substitutionQV,  len));
    EXPECT_TRUE(tests::CmpData<Nucleotide *>(seq.substitutionTag, seq2.substitutionTag, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.preBaseFrames,    seq2.preBaseFrames, len));
    EXPECT_TRUE(tests::CmpData<HalfWord *>(seq.widthInFrames,    seq2.widthInFrames, len));
};
