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
#ifndef BAM2BAX_CONVERTER_IMPL_HPP
#define BAM2BAX_CONVERTER_IMPL_HPP
#include <iostream>
#include "MetadataWriter.h"
#include "Bam2BaxInternal.h"
#include <pbbam/BamFile.h>
#include <pbbam/EntireFileQuery.h>


template<class T_HDFWRITER>
bool Bam2BaxConverter<T_HDFWRITER>::ConvertFile(void) {

    std::string outfn = settings_.outputBaxFilename;
    std::string infn = settings_.subreadsBamFilename;

    if (infn.empty()) infn = settings_.polymeraseBamFilename;

    PacBio::BAM::BamFile bamfile(infn);
    PacBio::BAM::BamHeader bamheader = bamfile.Header();

    if (bamheader.ReadGroups().size() != 1) {
        AddErrorMessage("Bam file must contain reads from exactly one SMRTCell.");
        return 0;
    }
    PacBio::BAM::ReadGroupInfo rg = bamheader.ReadGroups()[0];

    // Write metadata.xml to parent directory of Bax.h5.
    if (not settings_.outputMetadataFilename.empty())
        MetadataWriter metaWriter_(settings_.outputMetadataFilename, 
                                   rg, 
                                   settings_.outputAnalysisDirname);

    // Construct AcqParams
    AcqParams acqParams(Bam2BaxDefaults::Bax_ScanData_AduGain,
                        Bam2BaxDefaults::Bax_ScanData_CameraGain,
                        Bam2BaxDefaults::Bax_ScanData_CameraType,
                        Bam2BaxDefaults::Bax_ScanData_HotStartFrame,
                        Bam2BaxDefaults::Bax_ScanData_LaserOnFrame);

    // Construct scandata.
    ScanData scandata(acqParams);
    scandata.PlatformID(Sequel) // assume sequel movie 
            //.MovieName(rg.MovieName()) // rg.MovieName is not trust worthy due to an upstream bug
            .MovieName(settings_.movieName)
            .WhenStarted(rg.Date())
            .RunCode(Bam2BaxDefaults::Bax_ScanData_RunCode)  // bam does not contain RunCode
            .NumFrames(Bam2BaxDefaults::Bax_ScanData_NumFrames) // bam does not contain NumFrames
            //.FrameRate(strtof(rg.FrameRateHz().c_str(), NULL))
            .FrameRate(Bam2BaxDefaults::Bax_ScanData_FrameRate) // Ignore bam header FrameRate.
            .SequencingKit(rg.SequencingKit())
            .BindingKit(rg.BindingKit())
            .BaseMap(settings_.baseMap);

    // FIXME: pbbam needs to provide an API which returns BaseFeatures in read group
    std::vector<PacBio::BAM::BaseFeature> qvs = settings_.ignoreQV ? std::vector<PacBio::BAM::BaseFeature>({}) : internal::QVEnumsInFirstRecord(bamfile);

    // Regions attribute RegionTypes, which defines supported region types in ORDER.
    std::vector<RegionType> regionTypes = RegionTypeAdapter::ToRegionTypes(Bam2BaxDefaults::Bax_Regions_RegionTypes);
       
    if (not settings_.subreadsBamFilename.empty() and 
        not settings_.scrapsBamFilename.empty()) {

        T_HDFWRITER writer(outfn, 
                           scandata, 
                           rg.BasecallerVersion(), 
                           qvs,
                           Bam2BaxDefaults::Bax_Regions_RegionTypes);

        // Stich subreads and scraps in order to reconstruct polymerase reads.
        PacBio::BAM::VirtualPolymeraseReader reader(settings_.subreadsBamFilename,
                                                    settings_.scrapsBamFilename);
        while(reader.HasNext()) {
            // FIXME: pbbam should not crash when reading internal pulse features.
            const PacBio::BAM::VirtualPolymeraseBamRecord & record = reader.Next();
            SMRTSequence smrt;
            smrt.Copy(record, true);
            std::vector<RegionAnnotation> ras = RegionsAdapter::ToRegionAnnotations(record, regionTypes);
            if (not writer.WriteOneZmw(smrt, ras) or not writer.Errors().empty()) { break; }
            writer.Flush();
        }
        if (not settings_.ignoreQV) writer.WriteFakeDataSets();
        for (auto error: writer.Errors()) { AddErrorMessage(error); }
    } else if (not settings_.polymeraseBamFilename.empty()) {
        // Read polymerase reads from polymerase.bam directly.
        T_HDFWRITER writer(outfn, 
                           scandata, 
                           rg.BasecallerVersion(), 
                           qvs,
                           Bam2BaxDefaults::Bax_Regions_RegionTypes);

        PacBio::BAM::EntireFileQuery query(bamfile);
        for (auto record: query) {
            SMRTSequence smrt;
            smrt.Copy(record, true);
            RegionAnnotation ra(record.HoleNumber(), 
                                RegionTypeAdapter::ToRegionTypeIndex(PacBio::BAM::VirtualRegionType::HQREGION, regionTypes),
                                0, 0, 0);
            std::vector<RegionAnnotation> ras({ra});
            if (not writer.WriteOneZmw(smrt, ras) or not writer.Errors().empty()) { break; }
        }
        if (not settings_.ignoreQV) writer.WriteFakeDataSets();
        for (auto error: writer.Errors()) { AddErrorMessage(error); }
    }
            
    return errors_.empty();
}
#endif
