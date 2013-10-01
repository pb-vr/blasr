/*
 * =====================================================================================
 *
 *       Filename:  CCSIterator_gtest.cpp
 *
 *    Description:  Test common/files/CCSIterator.h
 *
 *        Version:  1.0
 *        Created:  11/29/2012 04:51:02 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */


#include "gtest/gtest.h"
#include "files/CCSIterator.h"
#include "datastructures/reads/RegionTable.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "../testdata.h"

using namespace std;

class CCSIteratorTestFixture: public testing::Test {
public:
    CCSIteratorTestFixture() {
    }

    void SetUp() {
        fileName = baxFile2; 
        reader = new HDFRegionTableReader();
        ccs = new CCSSequence();
        rgn = new RegionTable();

        int rev = reader->Initialize(fileName);
        EXPECT_TRUE(rev);
        reader->ReadTable(*rgn);
        reader->Close();

        rgn->SortTableByHoleNumber();
    }

    void TearDown() {
        if (reader) delete reader;
        if (ccs) delete ccs;
        if (rgn) delete rgn;
    }

    ~CCSIteratorTestFixture() {
    }

    string fileName;
    HDFRegionTableReader * reader;
    CCSSequence * ccs;
    RegionTable * rgn;
    CCSIterator it;
};


TEST_F(CCSIteratorTestFixture, Initialize) {
 
}
