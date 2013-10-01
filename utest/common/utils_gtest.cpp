/*
 * =====================================================================================
 *
 *       Filename:  utils_gtest.cpp
 *
 *    Description:  Test common/utils.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:22:58 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "utils.h"

TEST(UTILS, CrucialOpen) {
    ifstream ifs;
    string filename = "/nonexistingdir/nonexistingfile";
    //Expected behavour: exit with a message.
    EXPECT_EXIT( CrucialOpen(filename, ifs, std::ios::in),
        ::testing::ExitedWithCode(1), "");
}

TEST(UTILS, CeilOfFraction) {
    EXPECT_EQ(CeilOfFraction<int>(7, 100), 1);
    EXPECT_EQ(CeilOfFraction<int>(100, 7), 15);
    EXPECT_EQ(CeilOfFraction<int>(100, 99), 2);
    EXPECT_EQ(CeilOfFraction<int>(100, 5), 20);
}
