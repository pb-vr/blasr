/*
 * =====================================================================================
 *
 *       Filename:  TitleTable_gtest.cpp
 *
 *    Description:  Test common/datastructures/metagenome/TitleTable.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 03:34:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "datastructures/metagenome/TitleTable.h"
#include "../../testdata.h"

TEST(TitleTable, Read) {
    TitleTable tt;
    string fn = titleTable1;
    tt.Read(fn);
    EXPECT_STREQ(tt.table[0], "ref1 description1");
    EXPECT_STREQ(tt.table[1], "ref2 description2");
    EXPECT_EQ(tt.tableLength, 2);
    tt.Free();
}

