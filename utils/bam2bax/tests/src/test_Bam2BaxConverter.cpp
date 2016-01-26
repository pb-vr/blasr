// Copyright (c) 2014, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Yuan Li

#include "TestData.h"
#include "TestConstants.h"

#include "reads/RegionTable.hpp"
#include "Bam2BaxConverter.h"
#include <pbbam/virtual/VirtualRegion.h>
#include <string>
#include <gtest/gtest.h>

using namespace std;

std::vector<RegionType> defaultRegionTypes = RegionTable::DefaultRegionTypes();
std::vector<RegionType> definedRegionTypes = {Insert, HQRegion, Adapter, BarCode};

TEST(HDFBaxWriter, IsConvertibleVirtualRegionType)
{
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::ADAPTER, defaultRegionTypes));
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::HQREGION, defaultRegionTypes));
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::SUBREAD, defaultRegionTypes));
    EXPECT_FALSE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::BARCODE, defaultRegionTypes));
    EXPECT_FALSE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::LQREGION, defaultRegionTypes));

    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::ADAPTER, definedRegionTypes));
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::HQREGION, definedRegionTypes));
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::SUBREAD, definedRegionTypes));
    EXPECT_TRUE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::BARCODE, definedRegionTypes));
    EXPECT_FALSE(RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType::LQREGION, definedRegionTypes));
}

TEST(Bam2BaxConverter, CreateRegionAnnotation_DefaultRegionTypes)
{
    int zmw = 12134, beginPos = 0, endPos = 100, score = 770;

    EXPECT_EQ(defaultRegionTypes.size(), 3);
    EXPECT_EQ(defaultRegionTypes[0],  Adapter);
    EXPECT_EQ(defaultRegionTypes[1],  Insert);
    EXPECT_EQ(defaultRegionTypes[2],  HQRegion);

    PacBio::BAM::VirtualRegion vr(PacBio::BAM::VirtualRegionType::ADAPTER, beginPos, endPos, score);
    RegionAnnotation ra = RegionsAdapter::ToRegionAnnotation(zmw, vr, defaultRegionTypes);

    EXPECT_EQ(ra.GetHoleNumber(), zmw);
    EXPECT_EQ(defaultRegionTypes[ra.GetTypeIndex()], Adapter);
    EXPECT_EQ(ra.GetTypeString(defaultRegionTypes), "Adapter");
    EXPECT_EQ(ra.GetStart(),      beginPos);
    EXPECT_EQ(ra.GetEnd(),        endPos);
    EXPECT_EQ(ra.GetScore(),      score);

    PacBio::BAM::VirtualRegion vr2(PacBio::BAM::VirtualRegionType::HQREGION, beginPos, endPos, score);
    RegionAnnotation ra2 = RegionsAdapter::ToRegionAnnotation(zmw, vr2, defaultRegionTypes);
    EXPECT_EQ(defaultRegionTypes[ra2.GetTypeIndex()], HQRegion); 
    EXPECT_EQ(ra2.GetTypeString(defaultRegionTypes), "HQRegion");

    PacBio::BAM::VirtualRegion vr3(PacBio::BAM::VirtualRegionType::SUBREAD, beginPos, endPos, score);
    RegionAnnotation ra3 = RegionsAdapter::ToRegionAnnotation(zmw, vr3, defaultRegionTypes);
    EXPECT_EQ(defaultRegionTypes[ra3.GetTypeIndex()], Insert); 
    EXPECT_EQ(ra3.GetTypeString(defaultRegionTypes), "Insert");

};


TEST(Bam2BaxConverter, CreateRegionAnnotation_DefinedRegionTypes)
{
    int zmw = 12134, beginPos = 0, endPos = 100, score = 770;

    EXPECT_EQ(definedRegionTypes.size(), 4);
    EXPECT_EQ(definedRegionTypes[0],  Insert);
    EXPECT_EQ(definedRegionTypes[1],  HQRegion);
    EXPECT_EQ(definedRegionTypes[2],  Adapter);
    EXPECT_EQ(definedRegionTypes[3],  BarCode);

    PacBio::BAM::VirtualRegion vr(PacBio::BAM::VirtualRegionType::ADAPTER, beginPos, endPos, score);
    RegionAnnotation ra = RegionsAdapter::ToRegionAnnotation(zmw, vr, definedRegionTypes);

    EXPECT_EQ(ra.GetHoleNumber(), zmw);
    EXPECT_EQ(definedRegionTypes[ra.GetTypeIndex()], Adapter);
    EXPECT_EQ(ra.GetTypeString(definedRegionTypes), "Adapter");
    EXPECT_EQ(ra.GetStart(),      beginPos);
    EXPECT_EQ(ra.GetEnd(),        endPos);
    EXPECT_EQ(ra.GetScore(),      score);

    PacBio::BAM::VirtualRegion vr2(PacBio::BAM::VirtualRegionType::HQREGION, beginPos, endPos, score);
    RegionAnnotation ra2 = RegionsAdapter::ToRegionAnnotation(zmw, vr2, definedRegionTypes);
    EXPECT_EQ(definedRegionTypes[ra2.GetTypeIndex()], HQRegion); 
    EXPECT_EQ(ra2.GetTypeString(definedRegionTypes), "HQRegion");

    PacBio::BAM::VirtualRegion vr3(PacBio::BAM::VirtualRegionType::SUBREAD, beginPos, endPos, score);
    RegionAnnotation ra3 = RegionsAdapter::ToRegionAnnotation(zmw, vr3, definedRegionTypes);
    EXPECT_EQ(definedRegionTypes[ra3.GetTypeIndex()], Insert); 
    EXPECT_EQ(ra3.GetTypeString(definedRegionTypes), "Insert");

};
