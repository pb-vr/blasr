// Author: Yuan Li

#include "TestData.h"
#include "TestConstants.h"

#include "SMRTSequence.hpp"
#include "HDFGroup.hpp"
#include "HDFFile.hpp"
#include "HDFZMWMetricsWriter.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

TEST(HDFZMWMetricsWriter, EndToEnd)
{
    // setup
    std::string outfn = tests::Out_Dir  + "/" + "zmwmetrics.h5";

    HDFFile outfile;
    outfile.Open(outfn, H5F_ACC_TRUNC);

    HDFZMWMetricsWriter writer(outfn, outfile.rootGroup, tests::baseMap);

    SMRTSequence seq;
    tests::make_smrtseq(seq);

    bool OK = writer.WriteOneZmw(seq);

    EXPECT_TRUE(OK);
};

TEST(HDFZMWMetricsWriter, EndToEnd_RandomBaseMap)
{
    // setup
    std::string outfn = tests::Out_Dir  + "/" + "zmwmetrics2.h5";

    HDFFile outfile;
    outfile.Open(outfn, H5F_ACC_TRUNC);

    HDFZMWMetricsWriter writer(outfn, outfile.rootGroup, tests::randomBaseMap);

    SMRTSequence seq;
    tests::make_smrtseq(seq);

    bool OK = writer.WriteOneZmw(seq);

    EXPECT_TRUE(OK);
};
