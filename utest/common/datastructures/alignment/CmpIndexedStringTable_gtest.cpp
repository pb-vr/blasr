/*
 * =====================================================================================
 *
 *       Filename:  CmpIndexedStringTable_gtest.cpp
 *
 *    Description:  Test common/datastructures/alignment/CmpIndexedStringTable.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 01:50:24 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include "gtest/gtest.h"
#include "datastructures/alignment/CmpIndexedStringTable.h"

int    ids[5]   = {1,    2,    3,    30,    5};
string names[5] = {"n1", "n2", "n3", "n30", "n5"};
class CmpIndexedStringTableTest : public ::testing::Test {
public:
    //Be careful, SetUp() not Setup()
    virtual void SetUp() {
        for(int i = 0; i < 5; i++) {
            cmpIndexedStringTableTest.ids.push_back(ids[i]);
            cmpIndexedStringTableTest.names.push_back(names[i]);
        }
        cmpIndexedStringTableTest.StoreArrayIndexMap();

    }

    virtual void TearDown() {
    
    }

    CmpIndexedStringTable cmpIndexedStringTableTest;
};


//
// Test resize(int size)
//
TEST_F (CmpIndexedStringTableTest, resize) {
    cmpIndexedStringTableTest.resize(10);

    EXPECT_EQ(cmpIndexedStringTableTest.names.size(), 10);
    EXPECT_EQ(cmpIndexedStringTableTest.ids.size(), 10);
}

//
// Test StoreArrayIndexMap()
//
TEST_F (CmpIndexedStringTableTest, StoreArrayIndexMap) {
    for(int i = 0; i < cmpIndexedStringTableTest.ids.size(); i++) {
        EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(cmpIndexedStringTableTest.ids[i])->second, i);
    }
    EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(1 )->second,  0);
    EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(2 )->second,  1);
    EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(3 )->second,  2);
    EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(30)->second, 3);
    EXPECT_EQ(cmpIndexedStringTableTest.idToArrayIndex.find(5 )->second,  4);
}

//
// Test GetIndexOfId(int id, int & index)
//
TEST_F (CmpIndexedStringTableTest, GetIndexOfId) {
    int index = -1;
    int i = 0;
    for(i = 0; i < cmpIndexedStringTableTest.ids.size(); i++) {
        cmpIndexedStringTableTest.GetIndexOfId(ids[i], index);
        EXPECT_EQ(index, i);
    }
    i = 1;
    cmpIndexedStringTableTest.GetIndexOfId(i, index);
    EXPECT_EQ(index, 0);

    i = 2;
    cmpIndexedStringTableTest.GetIndexOfId(i, index);
    EXPECT_EQ(index, 1);

    i = 3;
    cmpIndexedStringTableTest.GetIndexOfId(i, index);
    EXPECT_EQ(index, 2);
    
    i = 5;
    cmpIndexedStringTableTest.GetIndexOfId(i, index);
    EXPECT_EQ(index, 4);
    
    i = 30;
    cmpIndexedStringTableTest.GetIndexOfId(i, index);
    EXPECT_EQ(index, 3);
}

//
// Test GetNameAtIndex(int index, string & name)
//
TEST_F (CmpIndexedStringTableTest, GetNameAtIndex) {
    // "Warning: the terminology of CmpIndexedStringTable.GetNameAtIndex
    // is confusing" 
        
    int index;
    string name;
    bool found;
    
    index = 1;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_EQ(name, "n1");
    EXPECT_TRUE(found);


    index = 2;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_EQ(name, "n2");
    EXPECT_TRUE(found);


    index = 3;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_EQ(name, "n3");
    EXPECT_TRUE(found);

    index = 30;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_EQ(name, "n30");
    EXPECT_TRUE(found);
    
    index = 5;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_EQ(name, "n5");
    EXPECT_TRUE(found);

    //
    //test getting name at an out-of-boundary index
    //
    index = 100;
    found = cmpIndexedStringTableTest.GetNameAtIndex(index, name);
    EXPECT_FALSE(found);
}

