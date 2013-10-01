/*
 * =====================================================================================
 *
 *       Filename:  defs_gtest.cpp
 *
 *    Description:  Test common/defs.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:17:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "defs.h"

TEST(DefsTest, MIN) {
    EXPECT_EQ(MIN(1,10000), 1);
    EXPECT_EQ(MIN(-1,10000), -1);
    EXPECT_EQ(MIN(-1,-2), -2);
}

TEST(DefsTest, MAX) {
    EXPECT_EQ(MAX(1,10000), 10000);
    EXPECT_EQ(MAX(-1,10000), 10000);
    EXPECT_EQ(MAX(-1,-2), -1);
}

TEST(DefsTest, SWAP) {
    int x = 10;
    int y = 100;
    SWAP(x, y);

    EXPECT_EQ(x, 100);
    EXPECT_EQ(y, 10);
}



