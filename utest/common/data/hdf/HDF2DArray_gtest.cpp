/*
 * ==========================================================================
 *
 *       Filename:  HDF2DArray_gtest.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  08/21/2013 07:00:49 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ==========================================================================
 */
#include <string>
#include "H5Cpp.h"
#include "data/hdf/HDFArray.h"
#include "data/hdf/HDFGroup.h"
#include "data/hdf/HDF2DArray.h"
#include "datastructures/reads/HoleXY.h"
#include "../../testdata.h"
#include "gtest/gtest.h"

using namespace std;
using namespace H5;
/*
//Note ::testing::Test not ::testing::TEST
//SetUp() and TearDown(), not Setup() and Teardown()
class HDF2DArrayTest: public :: testing::Test {
public:
    virtual void SetUp() { 
    }
    virtual void TearDown() {
    }
    string x;
};

TEST_F(HDF2DArrayTest, x) {
    EXPECT_TRUE(true);
}
*/

class HDF2DArrayTest : public ::testing::Test {
public:
    virtual void SetUp() {
        fileName = baxFile2;  // Defined in testdata.h
        try {
            FileAccPropList propList;
            pbihdfFile.openFile(fileName.c_str(), H5F_ACC_RDONLY, propList);
        } catch (Exception &e) {
            cout << "ERROR, could not open hdf file" << fileName
                << ", exiting." << endl;
            exit(1);
        }

        ASSERT_NE(rootGroup.Initialize(pbihdfFile, "/"), 0);
        ASSERT_NE(pulseDataGroup.Initialize(rootGroup, "PulseData"), 0);
        ASSERT_NE(baseCallsGroup.Initialize(pulseDataGroup, "BaseCalls"), 0);
        ASSERT_NE(zmwGroup.Initialize(baseCallsGroup, "ZMW"), 0);
    }

    virtual void TearDown() {
        rootGroup.Close();
        pulseDataGroup.Close();
        baseCallsGroup.Close();
        zmwGroup.Close();
        pbihdfFile.close();
    }

    string fileName;
    H5File pbihdfFile;
    HDFGroup rootGroup, pulseDataGroup, baseCallsGroup, zmwGroup;
};

//Test HDF2DArray, int16 
TEST_F(HDF2DArrayTest, int16) {
    HDF2DArray<int16_t> xyArray;
    xyArray.Initialize(zmwGroup, "HoleXY");
    vector<HoleXY> xys;
    xys.resize(10);
    for(int i = 0; i < 10; i++){
        xyArray.Read(i, i+1, xys[i].xy);
        ASSERT_EQ(xys[i].xy[0], 43);
        ASSERT_EQ(xys[i].xy[1], 44 + i);
    }

    int last = 54493;

    xyArray.Read(last, last+1, xys[0].xy);
    ASSERT_EQ(xys[0].xy[0], -43);
    ASSERT_EQ(xys[0].xy[1], -44);
    xyArray.Close();
}
