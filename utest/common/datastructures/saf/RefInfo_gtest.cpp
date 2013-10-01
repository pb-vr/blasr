/*
 * =====================================================================================
 *
 *       Filename:  RefInfo_gtest.cpp
 *
 *    Description:  Test common/datastructures/saf/RefInfo.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:02:00 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "datastructures/saf/RefInfo.h"

TEST(RefInfoTEST, RefIdToIndex) {
    RefInfo ri;

    string names[] = {"ref1", "ref3 description", "ref4 chr1 1:100"};
    ri.fullName.insert(ri.fullName.begin(), names, names+3);

    int id[] = {10, 7, 5};
    //int lens [] = {400, 300, 100};
    ri.id.insert(ri.id.begin(), id, id+3);

    int index;
    EXPECT_TRUE(ri.RefIdToIndex(10, index));
    EXPECT_EQ(index, 0);

    EXPECT_TRUE(ri.RefIdToIndex(7, index));
    EXPECT_EQ(index, 1);

    EXPECT_TRUE(ri.RefIdToIndex(5, index));
    EXPECT_EQ(index, 2);
}


