// Author: Yuan Li

#include "TestData.h"
#include "TestUtils.h"

#include "DNASequence.hpp"
#include <string>
#include <gtest/gtest.h>

using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

TEST(BAM2BAXTEST, EndToEnd)
{
    const std::string movieName = "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0";
    const std::string subreadsBamFilename = tests::Data_Dir + "/" + movieName + ".1.subreads.bam";
    const std::string scrapsBamFilename = tests::Data_Dir + "/" + movieName + ".1.scraps.bam";

    vector<string> bamFilenames = {subreadsBamFilename, scrapsBamFilename};
    const int result = RunBam2Bax(bamFilenames, "-o " + tests::Out_Dir + "/" + movieName + ".1");
    EXPECT_EQ(0, result);
}
