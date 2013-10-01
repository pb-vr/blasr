/*
 * =====================================================================================
 *
 *       Filename:  FileUtils_gtest.cpp
 *
 *    Description:  Test common/FileUtils.h
 *
 *        Version:  1.0
 *        Created:  10/29/2012 05:20:43 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */

#include "gtest/gtest.h"
#include "FileUtils.h"

string nonexistfile = "/nonexistingdir/nonexistingfile";
string readablefile = "/bin/ls";
string writeablefile = "/tmp/writabletmpfile";
string expected_errmsg = "Could not open file: " + nonexistfile;

TEST(FILEUTILS, CriticalOpenRead) {
    ifstream ifs;
    EXPECT_EXIT( CriticalOpenRead(nonexistfile, ifs, std::ios::in),
        ::testing::ExitedWithCode(1), expected_errmsg);
    CriticalOpenRead(readablefile, ifs, std::ios::in);
}

TEST(FILEUTILS, OpenRead) {
    ifstream ifs;
    EXPECT_EQ( OpenRead(nonexistfile, ifs, std::ios::in), 0);
    EXPECT_EQ( OpenRead(readablefile, ifs, std::ios::in), 1);
}


TEST(FILEUTILS, CriticalOpenWrite) {
    ofstream ofs;
    EXPECT_EXIT( CriticalOpenWrite(nonexistfile, ofs, std::ios::out),
        ::testing::ExitedWithCode(1), expected_errmsg);
    CriticalOpenWrite(writeablefile, ofs, std::ios::out);
}

TEST(FILEUTILS, OpenWrite) {
    ofstream ofs;
    EXPECT_EQ(OpenWrite(nonexistfile, ofs, std::ios::out), 0);
    EXPECT_EQ(OpenWrite(writeablefile, ofs, std::ios::out), 1);
}
