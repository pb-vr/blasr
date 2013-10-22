/*
 * =====================================================================================
 *
 *       Filename:  ScanData_gtest.cpp
 *
 *    Description:  Test common/datastructures/reads/ScanData.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 03:55:46 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "datastructures/reads/ScanData.h"

TEST(ScanDataTest, GetMovieName) {
    ScanData sd;
    EXPECT_EQ(sd.platformId, NoPlatform);
    EXPECT_EQ(sd.frameRate, 0);
    EXPECT_EQ(sd.numFrames, 0);
    EXPECT_EQ(sd.movieName.size(), 0);
    EXPECT_EQ(sd.runCode.size(), 0);
    EXPECT_EQ(sd.whenStarted.size(), 0);
    EXPECT_EQ(sd.baseMap.size(), 0);
}
