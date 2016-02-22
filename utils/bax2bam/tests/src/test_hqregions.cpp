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

// Author: Derek Barnett

#include "TestData.h"
#include "TestUtils.h"

#include "HDFBasReader.hpp"
#include "alignment/utils/RegionUtils.hpp"
#include "hdf/HDFRegionTableReader.hpp"
#include <boost/scoped_ptr.hpp>
#include <gtest/gtest.h>
#include <pbbam/BamFile.h>
#include <pbbam/BamRecord.h>
#include <pbbam/EntireFileQuery.h>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace PacBio;
using namespace PacBio::BAM;

TEST(HqRegionsTest, EndToEnd_Single)
{
    // setup
    const string movieName = "m140905_042212_sidney_c100564852550000001823085912221377_s1_X0";

    vector<string> baxFilenames;
    baxFilenames.push_back(tests::Data_Dir + "/" + movieName + ".1.bax.h5");

    const string generatedBam = movieName + ".hqregions.bam";
    const string scrapBam = movieName + ".lqregions.bam";

    // run conversion
    const int result = RunBax2Bam(baxFilenames, "--hqregion");
    EXPECT_EQ(0, result);

    {   // ensure PBIs exist
        const BamFile generatedBamFile(generatedBam);
        const BamFile scrapsBamFile(scrapBam);
        EXPECT_TRUE(generatedBamFile.PacBioIndexExists());
        EXPECT_TRUE(scrapsBamFile.PacBioIndexExists());
    }

    // open BAX reader on original data
    HDFBasReader baxReader;
    baxReader.IncludeField("Basecall");
    baxReader.IncludeField("DeletionQV");
    baxReader.IncludeField("DeletionTag");
    baxReader.IncludeField("InsertionQV");
    baxReader.IncludeField("PreBaseFrames");
    baxReader.IncludeField("MergeQV");
    baxReader.IncludeField("SubstitutionQV");
    baxReader.IncludeField("HQRegionSNR");
    // not using SubTag or PulseWidth here

    string baxBasecallerVersion;
    string baxBindingKit;
    string baxSequencingKit;

    const int initOk = baxReader.Initialize(baxFilenames.front());
    EXPECT_EQ(1, initOk);
    if (initOk == 1) {

        if (baxReader.scanDataReader.fileHasScanData && baxReader.scanDataReader.initializedRunInfoGroup) {

            if (baxReader.scanDataReader.runInfoGroup.ContainsAttribute("BindingKit")) {
                HDFAtom<std::string> bkAtom;
                if (bkAtom.Initialize(baxReader.scanDataReader.runInfoGroup, "BindingKit")) {
                    bkAtom.Read(baxBindingKit);
                    bkAtom.dataspace.close();
                }
            }

            if (baxReader.scanDataReader.runInfoGroup.ContainsAttribute("SequencingKit")) {
                HDFAtom<std::string> skAtom;
                if (skAtom.Initialize(baxReader.scanDataReader.runInfoGroup, "SequencingKit")) {
                    skAtom.Read(baxSequencingKit);
                    skAtom.dataspace.close();
                }
            }
        }

        baxReader.GetChangeListID(baxBasecallerVersion);
    }
        
    // compare 1st record from each file
    SMRTSequence baxRecord;
    EXPECT_TRUE(baxReader.GetNext(baxRecord) > 0);

    // read region table info
    boost::scoped_ptr<HDFRegionTableReader> regionTableReader(new HDFRegionTableReader);
    RegionTable regionTable;
    std::string fn = baxFilenames.front();
    EXPECT_TRUE(regionTableReader->Initialize(fn) != 0);
    regionTable.Reset();
    regionTableReader->ReadTable(regionTable);
    regionTableReader->Close();

    // Test primary, i.e. hqregions.bam, output
    EXPECT_NO_THROW(
    {
        // open BAM file
        BamFile bamFile(generatedBam);

        // check BAM header information
        const BamHeader& header = bamFile.Header();
        EXPECT_EQ(string("1.5"),     header.Version());
        EXPECT_EQ(string("unknown"), header.SortOrder());
        EXPECT_EQ(string("3.0.1"),   header.PacBioBamVersion());
        EXPECT_TRUE(header.Sequences().empty());
        EXPECT_TRUE(header.Comments().empty());
        ASSERT_FALSE(header.Programs().empty());

        const vector<string> readGroupIds = header.ReadGroupIds();
        ASSERT_FALSE(readGroupIds.empty());
        const ReadGroupInfo& rg = header.ReadGroup(readGroupIds.front());

        string rawId = movieName + "//HQREGION";
        string md5Id;
        MakeMD5(rawId, md5Id, 8);
        EXPECT_EQ(md5Id, rg.Id());

        EXPECT_EQ(string("PACBIO"), rg.Platform());
        EXPECT_EQ(movieName, rg.MovieName());

        EXPECT_TRUE(rg.SequencingCenter().empty());
        EXPECT_TRUE(rg.Date().empty());
        EXPECT_TRUE(rg.FlowOrder().empty());
        EXPECT_TRUE(rg.KeySequence().empty());
        EXPECT_TRUE(rg.Library().empty());
        EXPECT_TRUE(rg.Programs().empty());
        EXPECT_TRUE(rg.PredictedInsertSize().empty());
        EXPECT_TRUE(rg.Sample().empty());

        EXPECT_EQ("HQREGION", rg.ReadType());
        EXPECT_EQ(baxBasecallerVersion, rg.BasecallerVersion());
        EXPECT_EQ(baxBindingKit, rg.BindingKit());
        EXPECT_EQ(baxSequencingKit, rg.SequencingKit());
        EXPECT_EQ(75, std::stod(rg.FrameRateHz()));
        EXPECT_EQ("dq", rg.BaseFeatureTag(BaseFeature::DELETION_QV));
        EXPECT_EQ("dt", rg.BaseFeatureTag(BaseFeature::DELETION_TAG));
        EXPECT_EQ("iq", rg.BaseFeatureTag(BaseFeature::INSERTION_QV));
        EXPECT_EQ("ip", rg.BaseFeatureTag(BaseFeature::IPD));
        EXPECT_EQ("mq", rg.BaseFeatureTag(BaseFeature::MERGE_QV));
        EXPECT_EQ("sq", rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV));
        EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::SUBSTITUTION_TAG));
        EXPECT_EQ(FrameCodec::V1, rg.IpdCodec());

        int hqStart;
        int hqEnd;
        int hqScore;
        const bool lookupResult = LookupHQRegion(baxRecord.zmwData.holeNumber,
                                                 regionTable,
                                                 hqStart,
                                                 hqEnd,
                                                 hqScore);
        EXPECT_TRUE(lookupResult);

        vector<float> hqSnr;
        hqSnr.push_back(baxRecord.HQRegionSnr('A'));
        hqSnr.push_back(baxRecord.HQRegionSnr('C'));
        hqSnr.push_back(baxRecord.HQRegionSnr('G'));
        hqSnr.push_back(baxRecord.HQRegionSnr('T'));

        EXPECT_GT(hqSnr[0], 0);
        EXPECT_GT(hqSnr[1], 0);
        EXPECT_GT(hqSnr[2], 0);
        EXPECT_GT(hqSnr[3], 0);

        bool firstRecord = true;
        EntireFileQuery entireFile(bamFile);
        for ( BamRecord& bamRecord : entireFile) {
            if (!firstRecord)
                break;
            firstRecord = false;

            const BamRecordImpl& bamRecordImpl = bamRecord.Impl();
            EXPECT_EQ(4680,bamRecordImpl.Bin());
            EXPECT_EQ(0,   bamRecordImpl.InsertSize());
            EXPECT_EQ(255, bamRecordImpl.MapQuality());
            EXPECT_EQ(-1,  bamRecordImpl.MatePosition());
            EXPECT_EQ(-1,  bamRecordImpl.MateReferenceId());
            EXPECT_EQ(-1,  bamRecordImpl.Position());
            EXPECT_EQ(-1,  bamRecordImpl.ReferenceId());
            EXPECT_FALSE(bamRecordImpl.IsMapped());

            const int holeNumber    = baxRecord.zmwData.holeNumber;
            const int subreadStart  = hqStart;
            const int subreadEnd    = hqEnd;

            const string expectedName = movieName + "/" +
                    to_string(holeNumber)   + "/" +
                    to_string(subreadStart) + "_" +
                    to_string(subreadEnd);
            EXPECT_EQ(expectedName, bamRecordImpl.Name());

            using PacBio::BAM::QualityValue;
            using PacBio::BAM::QualityValues;

            const DNALength length = subreadEnd - subreadStart;

            string expectedSequence;
            expectedSequence.assign((const char*)baxRecord.seq + subreadStart, length);

            const string bamSequence = bamRecord.Sequence();
            const QualityValues bamQualities = bamRecord.Qualities();
            EXPECT_EQ(expectedSequence, bamSequence);
            EXPECT_TRUE(bamQualities.empty());

            const QualityValues bamDeletionQVs = bamRecord.DeletionQV();
            const QualityValues bamInsertionQVs = bamRecord.InsertionQV();
            const QualityValues bamMergeQVs = bamRecord.MergeQV();
            const QualityValues bamSubstitutionQVs = bamRecord.SubstitutionQV();

            for (size_t i = 0; i < length; ++i) {
                const size_t pos = subreadStart + i;

                EXPECT_EQ((QualityValue)baxRecord.GetDeletionQV(pos),     bamDeletionQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetInsertionQV(pos),    bamInsertionQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetMergeQV(pos),        bamMergeQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetSubstitutionQV(pos), bamSubstitutionQVs.at(i));
            }

            if (baxRecord.deletionTag)
            {
                string expectedDeletionTags;
                expectedDeletionTags.assign((char*)baxRecord.deletionTag + subreadStart,
                                            (char*)baxRecord.deletionTag + subreadStart + length);
                const string& bamDeletionTags = bamRecord.DeletionTag();
                EXPECT_EQ(expectedDeletionTags, bamDeletionTags);
            }

            if (baxRecord.substitutionTag)
            {
                string expectedSubstitutionTags;
                expectedSubstitutionTags.assign((char*)baxRecord.substitutionTag + subreadStart,
                                            (char*)baxRecord.substitutionTag + subreadStart + length);
                const string& bamSubstitutionTags = bamRecord.SubstitutionTag();
                EXPECT_EQ(expectedSubstitutionTags, bamSubstitutionTags);
            }

            // TODO: IPDs

            EXPECT_EQ(md5Id,        bamRecord.ReadGroupId());
            EXPECT_EQ(movieName,    bamRecord.MovieName());
            EXPECT_EQ(1,            bamRecord.NumPasses());
            EXPECT_EQ(holeNumber,   bamRecord.HoleNumber());
            EXPECT_EQ(subreadStart, bamRecord.QueryStart());
            EXPECT_EQ(subreadEnd,   bamRecord.QueryEnd());
            EXPECT_EQ(hqSnr,        bamRecord.SignalToNoise());
            EXPECT_FALSE(bamRecord.HasLocalContextFlags());
        }

    }); // EXPECT_NO_THROW

    // Test secondary, i.e. lqregions.bam, output
    EXPECT_NO_THROW(
    {
        // open BAM file
        BamFile bamFile(scrapBam);

        // check BAM header information
        const BamHeader& header = bamFile.Header();
        EXPECT_EQ(string("1.5"),     header.Version());
        EXPECT_EQ(string("unknown"), header.SortOrder());
        EXPECT_EQ(string("3.0.1"),   header.PacBioBamVersion());
        EXPECT_TRUE(header.Sequences().empty());
        EXPECT_TRUE(header.Comments().empty());
        ASSERT_FALSE(header.Programs().empty());

        const vector<string> readGroupIds = header.ReadGroupIds();
        ASSERT_FALSE(readGroupIds.empty());
        const ReadGroupInfo& rg = header.ReadGroup(readGroupIds.front());

        string rawId = movieName + "//SCRAP";
        string md5Id;
        MakeMD5(rawId, md5Id, 8);
        EXPECT_EQ(md5Id, rg.Id());

        EXPECT_EQ(string("PACBIO"), rg.Platform());
        EXPECT_EQ(movieName, rg.MovieName());

        EXPECT_TRUE(rg.SequencingCenter().empty());
        EXPECT_TRUE(rg.Date().empty());
        EXPECT_TRUE(rg.FlowOrder().empty());
        EXPECT_TRUE(rg.KeySequence().empty());
        EXPECT_TRUE(rg.Library().empty());
        EXPECT_TRUE(rg.Programs().empty());
        EXPECT_TRUE(rg.PredictedInsertSize().empty());
        EXPECT_TRUE(rg.Sample().empty());

        EXPECT_EQ("SCRAP", rg.ReadType());
        EXPECT_EQ(baxBasecallerVersion, rg.BasecallerVersion());
        EXPECT_EQ(baxBindingKit, rg.BindingKit());
        EXPECT_EQ(baxSequencingKit, rg.SequencingKit());
        EXPECT_EQ(75, std::stod(rg.FrameRateHz()));
        EXPECT_EQ("dq", rg.BaseFeatureTag(BaseFeature::DELETION_QV));
        EXPECT_EQ("dt", rg.BaseFeatureTag(BaseFeature::DELETION_TAG));
        EXPECT_EQ("iq", rg.BaseFeatureTag(BaseFeature::INSERTION_QV));
        EXPECT_EQ("ip", rg.BaseFeatureTag(BaseFeature::IPD));
        EXPECT_EQ("mq", rg.BaseFeatureTag(BaseFeature::MERGE_QV));
        EXPECT_EQ("sq", rg.BaseFeatureTag(BaseFeature::SUBSTITUTION_QV));
        EXPECT_FALSE(rg.HasBaseFeature(BaseFeature::SUBSTITUTION_TAG));
        EXPECT_EQ(FrameCodec::V1, rg.IpdCodec());

        int hqStart;
        int hqEnd;
        int hqScore;
        const bool lookupResult = LookupHQRegion(baxRecord.zmwData.holeNumber,
                                                 regionTable,
                                                 hqStart,
                                                 hqEnd,
                                                 hqScore);
        EXPECT_TRUE(lookupResult);

        vector<float> hqSnr;
        hqSnr.push_back(baxRecord.HQRegionSnr('A'));
        hqSnr.push_back(baxRecord.HQRegionSnr('C'));
        hqSnr.push_back(baxRecord.HQRegionSnr('G'));
        hqSnr.push_back(baxRecord.HQRegionSnr('T'));

        EXPECT_GT(hqSnr[0], 0);
        EXPECT_GT(hqSnr[1], 0);
        EXPECT_GT(hqSnr[2], 0);
        EXPECT_GT(hqSnr[3], 0);

        bool firstRecord = true;
        EntireFileQuery entireFile(bamFile);
        for ( BamRecord& bamRecord : entireFile) {
            if (!firstRecord)
                break;
            firstRecord = false;

            const BamRecordImpl& bamRecordImpl = bamRecord.Impl();
            EXPECT_EQ(4680,bamRecordImpl.Bin());
            EXPECT_EQ(0,   bamRecordImpl.InsertSize());
            EXPECT_EQ(255, bamRecordImpl.MapQuality());
            EXPECT_EQ(-1,  bamRecordImpl.MatePosition());
            EXPECT_EQ(-1,  bamRecordImpl.MateReferenceId());
            EXPECT_EQ(-1,  bamRecordImpl.Position());
            EXPECT_EQ(-1,  bamRecordImpl.ReferenceId());
            EXPECT_FALSE(bamRecordImpl.IsMapped());

            const int holeNumber    = baxRecord.zmwData.holeNumber;
            const int subreadStart  = 0;
            const int subreadEnd    = hqStart;

            const string expectedName = movieName + "/" +
                    to_string(holeNumber)   + "/" +
                    to_string(subreadStart) + "_" +
                    to_string(subreadEnd);
            EXPECT_EQ(expectedName, bamRecordImpl.Name());

            using PacBio::BAM::QualityValue;
            using PacBio::BAM::QualityValues;

            const DNALength length = subreadEnd - subreadStart;

            string expectedSequence;
            expectedSequence.assign((const char*)baxRecord.seq + subreadStart, length);

            const string bamSequence = bamRecord.Sequence();
            const QualityValues bamQualities = bamRecord.Qualities();
            EXPECT_EQ(expectedSequence, bamSequence);
            EXPECT_TRUE(bamQualities.empty());

            const QualityValues bamDeletionQVs = bamRecord.DeletionQV();
            const QualityValues bamInsertionQVs = bamRecord.InsertionQV();
            const QualityValues bamMergeQVs = bamRecord.MergeQV();
            const QualityValues bamSubstitutionQVs = bamRecord.SubstitutionQV();

            for (size_t i = 0; i < length; ++i) {
                const size_t pos = subreadStart + i;

                EXPECT_EQ((QualityValue)baxRecord.GetDeletionQV(pos),     bamDeletionQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetInsertionQV(pos),    bamInsertionQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetMergeQV(pos),        bamMergeQVs.at(i));
                EXPECT_EQ((QualityValue)baxRecord.GetSubstitutionQV(pos), bamSubstitutionQVs.at(i));
            }

            if (baxRecord.deletionTag)
            {
                string expectedDeletionTags;
                expectedDeletionTags.assign((char*)baxRecord.deletionTag + subreadStart,
                                            (char*)baxRecord.deletionTag + subreadStart + length);
                const string& bamDeletionTags = bamRecord.DeletionTag();
                EXPECT_EQ(expectedDeletionTags, bamDeletionTags);
            }

            if (baxRecord.substitutionTag)
            {
                string expectedSubstitutionTags;
                expectedSubstitutionTags.assign((char*)baxRecord.substitutionTag + subreadStart,
                                            (char*)baxRecord.substitutionTag + subreadStart + length);
                const string& bamSubstitutionTags = bamRecord.SubstitutionTag();
                EXPECT_EQ(expectedSubstitutionTags, bamSubstitutionTags);
            }

            // TODO: IPDs
            EXPECT_EQ(md5Id,        bamRecord.ReadGroupId());
            EXPECT_EQ(movieName,    bamRecord.MovieName());
            EXPECT_EQ(1,            bamRecord.NumPasses());
            EXPECT_EQ(holeNumber,   bamRecord.HoleNumber());
            EXPECT_EQ(subreadStart, bamRecord.QueryStart());
            EXPECT_EQ(subreadEnd,   bamRecord.QueryEnd());
            EXPECT_EQ(hqSnr,        bamRecord.SignalToNoise());
            EXPECT_FALSE(bamRecord.HasLocalContextFlags());
        }

    }); // EXPECT_NO_THROW

    // cleanup
    EXPECT_NO_THROW(
    {
        baxReader.Close();
        RemoveFile(generatedBam);
        RemoveFile(scrapBam);
        RemoveFile(generatedBam + ".pbi");
        RemoveFile(scrapBam + ".pbi");
    });
}
