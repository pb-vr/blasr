// Copyright (c) 2014-2015, Pacific Biosciences of California, Inc.
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

#ifndef _BAM2BAXINTERNAL_H_
#define _BAM2BAXINTERNAL_H_

//namespace internal
namespace internal {
    /// \name \{
    static const std::vector<PacBio::BAM::BaseFeature> QVEnums = {
          PacBio::BAM::BaseFeature::DELETION_QV
        , PacBio::BAM::BaseFeature::DELETION_TAG
        , PacBio::BAM::BaseFeature::INSERTION_QV
        , PacBio::BAM::BaseFeature::MERGE_QV
        , PacBio::BAM::BaseFeature::SUBSTITUTION_QV
        , PacBio::BAM::BaseFeature::SUBSTITUTION_TAG
        , PacBio::BAM::BaseFeature::IPD
        , PacBio::BAM::BaseFeature::PULSE_WIDTH
        , PacBio::BAM::BaseFeature::PKMID
        , PacBio::BAM::BaseFeature::PKMEAN
        , PacBio::BAM::BaseFeature::LABEL
        , PacBio::BAM::BaseFeature::LABEL_QV
        , PacBio::BAM::BaseFeature::ALT_LABEL
        , PacBio::BAM::BaseFeature::ALT_LABEL_QV
        , PacBio::BAM::BaseFeature::PULSE_MERGE_QV
        , PacBio::BAM::BaseFeature::PULSE_CALL
        , PacBio::BAM::BaseFeature::START_FRAME
        , PacBio::BAM::BaseFeature::PULSE_CALL_WIDTH
    };

    /// \returns QVs contained by read group rg.
    /// FIXME: this function should be provided by pbbam.ReadGroupInfo
    /// FIXME: pbbam, ReadGroupInfo does not recognize internal pulse features such as AltLabelQV.
    inline std::vector<PacBio::BAM::BaseFeature> 
    QVEnumsInReadGroup(const PacBio::BAM::ReadGroupInfo & rg) {
        std::vector<PacBio::BAM::BaseFeature> ret;
        for (auto it = internal::QVEnums.begin(); it != internal::QVEnums.end(); it++) {
            if (rg.HasBaseFeature(*it)) {
                ret.push_back(*it);
            }
        }
        return ret;
    }
    /// \}
 
    /// \returns QVs contained by the first record if it exists, otherwise, return {}
    /// FIXME: this function provides an alternative route to get QVs contained in the bam file now,
    /// because pbbam ReadGroupInfo does not recorgize internal pulse features such as AltLabelQV.
    /// Note: Ignore Label because it is neither base feature nor internal pulse feature.
    inline std::vector<PacBio::BAM::BaseFeature> 
    QVEnumsInFirstRecord(const PacBio::BAM::BamFile & bamFile) {
        std::vector<PacBio::BAM::BaseFeature> ret;
        PacBio::BAM::EntireFileQuery query(bamFile);
        for (const PacBio::BAM::BamRecord & record: query) {
            if (record.HasDeletionQV())      {ret.push_back(PacBio::BAM::BaseFeature::DELETION_QV);}
            if (record.HasDeletionTag())     {ret.push_back(PacBio::BAM::BaseFeature::DELETION_TAG);}
            if (record.HasInsertionQV())     {ret.push_back(PacBio::BAM::BaseFeature::INSERTION_QV);}
            if (record.HasMergeQV())         {ret.push_back(PacBio::BAM::BaseFeature::MERGE_QV);}
            if (record.HasSubstitutionQV())  {ret.push_back(PacBio::BAM::BaseFeature::SUBSTITUTION_QV);}
            if (record.HasSubstitutionTag()) {ret.push_back(PacBio::BAM::BaseFeature::SUBSTITUTION_TAG);}
            if (record.HasIPD())             {ret.push_back(PacBio::BAM::BaseFeature::IPD);}
            if (record.HasPulseWidth())      {ret.push_back(PacBio::BAM::BaseFeature::PULSE_WIDTH);}
            if (record.HasPkmid())           {ret.push_back(PacBio::BAM::BaseFeature::PKMID);}
            if (record.HasPkmean())          {ret.push_back(PacBio::BAM::BaseFeature::PKMEAN);}
            if (record.HasLabelQV())         {ret.push_back(PacBio::BAM::BaseFeature::LABEL_QV);}
            if (record.HasAltLabelTag())     {ret.push_back(PacBio::BAM::BaseFeature::ALT_LABEL);}
            if (record.HasAltLabelQV())      {ret.push_back(PacBio::BAM::BaseFeature::ALT_LABEL_QV);}
            if (record.HasPulseMergeQV())    {ret.push_back(PacBio::BAM::BaseFeature::PULSE_MERGE_QV);}
            if (record.HasPulseCall())       {ret.push_back(PacBio::BAM::BaseFeature::PULSE_CALL);}
            if (record.HasStartFrame())      {ret.push_back(PacBio::BAM::BaseFeature::START_FRAME);}
            if (record.HasPulseCallWidth())  {ret.push_back(PacBio::BAM::BaseFeature::PULSE_CALL_WIDTH);}
            break; // only use the first record.
        }
        return ret;
    }
};

#endif
