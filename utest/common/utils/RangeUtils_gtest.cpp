/*
 * =====================================================================================
 *
 *       Filename:  RangeUtils_gtest.cpp
 *
 *    Description:  Test common/utils/RangeUtils.h
 *
 *        Version:  1.0
 *        Created:  05/2/2013 06:01:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "utils/RangeUtils.h"

TEST(RangeTest, RangeConstructor) {
    Range r(1, 2);
    EXPECT_EQ(r.start, 1);
    EXPECT_EQ(r.end, 2);

    Range r2(1);
    EXPECT_EQ(r2.end, 1);
}

TEST(RangeTest, Ranges) {
    UInt queryInRange[11] = {1,2,3,4,10,11,12,13,14,15,20};
    UInt queryNotInRange[8] = {0, 16, 17, 18, 19, 30, 5, 100000}; 

    Ranges ranges1(string("1,2,3,4,10-15,20-20"));

    for (int i = 0; i < 11; i++) {
        EXPECT_TRUE(ranges1.contains(queryInRange[i]));
    }

    for (int i = 0; i < 8; i++) {
        EXPECT_FALSE(ranges1.contains(queryNotInRange[i]));
    }
}

TEST(RangeTest, SetRanges) {
    Ranges rs;
    rs.setRanges("199");

    EXPECT_TRUE(rs.contains(199));
    EXPECT_EQ(rs.size(), 1);
}


TEST(RangeTest, max) {
    Ranges rs("199");
    EXPECT_EQ(rs.max(), 199);
    Ranges rs1("1, 10000, 10-30, 4000-5000");
    EXPECT_EQ(rs1.max(), 10000);

    Ranges rs2("1, 1000, 10-30, 4000-5000, 633-877");
    EXPECT_EQ(rs2.max(), 5000);
}
