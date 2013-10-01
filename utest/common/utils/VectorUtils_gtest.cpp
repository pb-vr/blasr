/*
 * =====================================================================================
 *
 *       Filename:  VectorUtils_gtest.cpp
 *
 *    Description:  Test common/utils/VectorUtils.h
 *
 *        Version:  1.0
 *        Created:  01/17/2013 06:01:01 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "utils/VectorUtils.h"

// Test ClearMemory(vector<T> vt)
TEST(VectorUtils, ClearMemory) {
    vector<int> vi;
    vi.push_back(1);

    int size = 1000000;
    vi.reserve(size);
    EXPECT_EQ(vi.size(), 1);
    EXPECT_EQ(vi.capacity(), size);

    ClearMemory(vi);
    EXPECT_EQ(vi.size(), 0);
    EXPECT_EQ(vi.capacity(), 0);
} 

