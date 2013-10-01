/*
 * =====================================================================================
 *
 *       Filename:  RegionUtils_gtest.cpp
 *
 *    Description:  Test common/utils/RegionUtils.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 05:11:56 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "utils/RegionUtils.h"
#include "datastructures/reads/RegionTable.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "../testdata.h"

class RegionUtilTestFixture: public testing::Test {
public:
    RegionUtilTestFixture() {
        fileName = baxFile1; 
        reader = new HDFRegionTableReader();
        int rev = reader->Initialize(fileName);
        EXPECT_TRUE(rev);
        reader->ReadTable(regionTable);
        reader->Close();
    }

    void SetUp() {
        /* The region table for hole number 14798 is as follows.
         * type start end  score
         *   1   0    712   -1 
         *   1   760  2040  -1 
         *   1   2098 3452  -1
         *   0   712  760   937
         *   0   2040 2098  741
         *   2   0    3424  819
         *   where type 1 = insertion, 0 = adapter, 2 = HQRegion
         * */
        hqStart = 0;
        hqEnd   = 3424;
        hqScore = 819;
        holeNumber = 14798;
    }

    void TearDown() {
    }

    ~RegionUtilTestFixture() {
        delete reader;
    }

    HDFRegionTableReader * reader;
    RegionTable regionTable;
    string fileName;
    int hqStart, hqEnd, hqScore, holeNumber;
};

TEST_F(RegionUtilTestFixture, LookupHQRegion) {
    int start, end, score;
    bool rev = LookupHQRegion(holeNumber, regionTable, start, end, score);
    EXPECT_EQ(rev, true);
    EXPECT_EQ(start, hqStart);
    EXPECT_EQ(end, hqEnd);
    EXPECT_EQ(score, hqScore);
}

TEST_F(RegionUtilTestFixture, GetHighQulitySubreadsIntervals) {
    vector<ReadInterval> intervals;
    intervals.push_back(ReadInterval(0, 712));
    intervals.push_back(ReadInterval(760, 2040));
    intervals.push_back(ReadInterval(2098, 3452));

    vector<int> directions;
    directions.push_back(0);
    directions.push_back(1);
    directions.push_back(0);

    int indx = GetHighQualitySubreadsIntervals(intervals, directions, hqStart, hqEnd);
    EXPECT_EQ(intervals.size(), 3);
    EXPECT_EQ(indx, 2);
    int starts [3] = {0, 760, 2098};
    int ends   [3] = {712, 2040, 3424};
    int ds     [3] = {0, 1, 0};
    for(int i=0; i < 3; i++) {
        EXPECT_EQ(intervals[i].start, starts[i]);
        EXPECT_EQ(intervals[i].end  , ends[i]  );
        EXPECT_EQ(directions[i]     , ds[i]    );
    }

    indx = GetHighQualitySubreadsIntervals(intervals, directions, hqStart, hqEnd, 800);
    EXPECT_EQ(intervals.size(), 2);
    // The first interval and its direction has been removed as the length is less
    // than 800.
    for(int i=0; i < 2; i++) {
        EXPECT_EQ(intervals[i].start, starts[i+1]);
        EXPECT_EQ(intervals[i].end  , ends[i+1]  );
    }


}
