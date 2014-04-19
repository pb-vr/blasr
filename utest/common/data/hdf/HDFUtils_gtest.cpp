/*
 * ============================================================================
 *
 *       Filename:  HDFUtils_gtest.cpp
 *
 *    Description:  Test common/data/hdf/HDFUtils.h
 *
 *        Version:  1.0
 *        Created:  11/04/2013 10:55:10 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */

#include "data/hdf/HDFUtils.h"
#include "gtest/gtest.h"
#include "../../testdata.h"

TEST(HDFUtils, Create) {
    EXPECT_EQ(GetH5MovieName(baxFile1), 
        "m130220_114643_42129_c100471902550000001823071906131347_s1_p0");

    EXPECT_EQ(GetH5MovieName(baxFile2),
        "m130220_114643_42129_c100471902550000001823071906131347_s1_p0");

    EXPECT_EQ(GetH5MovieName(plsFile1),
        "m121215_065521_richard_c100425710150000001823055001121371_s1_p0");

    EXPECT_EQ(GetH5MovieName(ccsFile1),
        "m130328_211423_ethan_c100499512550000001823070408081371_s1_p0");

    EXPECT_EQ(GetH5MovieName(rgnFile1),
        "m130427_152935_42178_c100518602550000001823079209281316_s1_p0");
}
