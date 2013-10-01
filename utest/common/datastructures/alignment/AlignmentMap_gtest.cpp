/*
 * =====================================================================================
 *
 *       Filename:  AlignmentMap_gtest.cpp
 *
 *    Description:  Test common/datastructures/alignment/AlignmentMap.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 01:47:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "datastructures/alignment/AlignmentMap.h"
#include <climits>
#include <iostream>
#include <fstream>

// Test 
// void CreateSequenceToAlignmentMap(const string & alignedSequence,
//         vector<int> & baseToAlignmentMap); 
TEST(AlignmentMap, CreateSequenceToAlignmentMap) {
    const string & alignedSequence1 = "ATCTGAG-AAA-";
    const int size1 = 10;
    int map1[size1] = {0, 1, 2, 3, 4,5, 6, 8, 9, 10};

    const string & alignedSequence2 = "--ATCTGAG----AAA----";
    const int size2 = 10;
    int map2[size2] = {2,3,4,5,6,7,8,13,14,15};

    const string & alignedSequence3 = "-------";
    const int size3 = 0;

    vector <int> resMap1;
    CreateSequenceToAlignmentMap(alignedSequence1, resMap1); 
    EXPECT_EQ(size1, resMap1.size());
    for (int i = 0; i < size1; i++) {
        // cout << resMap1[i] << ", ";
        EXPECT_EQ(map1[i], resMap1[i]);
    }
    // cout << endl;

    vector <int> resMap2;
    CreateSequenceToAlignmentMap(alignedSequence2, resMap2);
    EXPECT_EQ(size2, resMap2.size());
    for (int i = 0; i < size2; i++) {
    //    cout << resMap1[i] << ", ";
        EXPECT_EQ(map2[i], resMap2[i]);
    }
    // cout << endl;
    
    vector <int> resMap3;
    CreateSequenceToAlignmentMap(alignedSequence3, resMap3);
    EXPECT_EQ(size3, resMap3.size());

}
