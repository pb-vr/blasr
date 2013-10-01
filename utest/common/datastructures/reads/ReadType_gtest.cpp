/*
 * ==================================================================
 *
 *       Filename:  ReadType_gtest.cpp
 *
 *    Description:  Test common/datastructures/reads/ReadType.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 03:54:55 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ==================================================================
 */
#include "gtest/gtest.h"
#include "datastructures/reads/ReadType.h"


TEST(ReadTypeTest, ParseReadType) {

    string standard   = "Standard";
    string ccs        = "CCS";
    string rccs       = "RCCS";
    string noreadtype = "standard";

    EXPECT_EQ(ParseReadType(standard), ReadType::Standard);
    EXPECT_EQ(ParseReadType(ccs),      ReadType::CCS);
    EXPECT_EQ(ParseReadType(rccs),     ReadType::RCCS);
    EXPECT_EQ(ParseReadType(noreadtype), ReadType::NoReadType);
}
