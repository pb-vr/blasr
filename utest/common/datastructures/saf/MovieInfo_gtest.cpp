/*
 * =====================================================================================
 *
 *       Filename:  MovieInfo_gtest.cpp
 *
 *    Description:  Test common/datastructures/saf/MovieInfo.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:01:21 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "Types.h"
#include "datastructures/saf/MovieInfo.h"

// Test MovieInfo.FindMovie
TEST(MovieInfoTest, FindMovie) {
    MovieInfo movieInfo;
    int ids[5] = {3, 5, 7, 4, 0} ;
    string names[5] = {"movieX", "movieY", "abc", "m000000000032102389170_s0", "m000000000032102389170_s"};

    for (int i = 0; i < 5; i++) {
        movieInfo.name.push_back(names[i]);
        movieInfo.id.push_back(ids[i]);
    }

    string nameVal, nameVal2;
    int ret = movieInfo.FindMovie(4, nameVal);
    EXPECT_EQ(nameVal, names[3]);
    EXPECT_EQ(ret, 1);

    ret = movieInfo.FindMovie(-1, nameVal);
    EXPECT_EQ(nameVal2, "");
    EXPECT_EQ(ret, 0);

}
