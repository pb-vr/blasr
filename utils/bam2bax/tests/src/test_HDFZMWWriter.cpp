// Author: Yuan Li

#include "TestData.h"

#include "SMRTSequence.hpp"
#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFZMWWriter.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

TEST(HDFZMWWriter, EndToEnd)
{
    // setup
    std::string outfn = tests::Out_Dir  + "/" + "zmw.h5";

    HDFFile outfile;
    outfile.Open(outfn, H5F_ACC_TRUNC);

    HDFZMWWriter writer(outfn, outfile.rootGroup);

    for (int i = 1 ; i < 1000; i++) {
        SMRTSequence seq;
        seq.length = i;
        seq.zmwData.holeNumber = i;
        seq.zmwData.holeStatus = static_cast<unsigned char> (8);
        seq.readScore = 0.87;
        seq.HQRegionSnr('A', 0.1);
        seq.HQRegionSnr('C', 0.2);
        seq.HQRegionSnr('G', 0.3);
        seq.HQRegionSnr('T', 0.4);

        bool OK = writer.WriteOneZmw(seq);
        EXPECT_TRUE(OK);
    }
};
