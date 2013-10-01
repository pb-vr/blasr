/*
 * =====================================================================================
 *
 *       Filename:  AlnGroup_gtest.cpp
 *
 *    Description:  Test common/datastructures/saf/AlnGroup.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:00:29 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include "gtest/gtest.h"
#include "datastructures/saf/AlnGroup.h"

//Test AlnGroup.FindPath();
TEST(AlnGroupTest, FindPath) {
    AlnGroup alnGroup;
    unsigned int ids[10] = {3, 4, 6, 1, 2, 8, 9, 12, 11, 7};
    string paths[10] = {"path1", "path2", "path3", "path4", "path5",
                        "path6", "path7", "path8", "path9", "path10"};
    for(int i = 0; i < 10; i++) {
        alnGroup.id.push_back(ids[i]);
        alnGroup.path.push_back(paths[i]);
    }

    string val, val1;
    int ret = alnGroup.FindPath(3, val);
    EXPECT_EQ(val, paths[0]);
    EXPECT_EQ(ret, 1);

    
    ret = alnGroup.FindPath(100, val1);
    EXPECT_EQ(ret, 0);
    EXPECT_EQ(val1, "");

}
