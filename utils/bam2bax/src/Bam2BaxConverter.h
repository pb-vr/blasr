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
#ifndef _BAM2BAXCONVERTER_H_
#define _BAM2BAXCONVERTER_H_

#include <string>
#include <vector>
#include <algorithm>
#include <pbbam/BamFile.h>
#include <pbbam/BamHeader.h>
#include <pbbam/ReadGroupInfo.h>
#include <pbbam/virtual/VirtualPolymeraseReader.h>
#include <pbbam/virtual/VirtualPolymeraseBamRecord.h>
#include <pbbam/virtual/VirtualRegion.h>
#include <pbbam/virtual/VirtualRegionType.h>
#include <pbbam/virtual/VirtualRegionTypeMap.h>
#include "RegionsAdapter.h"
#include "IConverter.h"
#include <boost/scoped_ptr.hpp>


namespace Bam2BaxDefaults {
    // Default value of attribute /ScanData/AcqParams/NumFrames in Bax.
    static const unsigned int Bax_ScanData_NumFrames = 0;
    // Default value of attribute /ScanData/AcqParams/AduGain in Bax.
    static const float Bax_ScanData_AduGain = 1.0;
    // Default value of attribute /ScanData/AcqParams/CameraGain in Bax.
    static const float Bax_ScanData_CameraGain = 1.0;
    // Default value of attribute /ScanData/AcqParams/CameraType in Bax.
    static const int Bax_ScanData_CameraType = 0;
    // Default value of attribute /ScanData/AcqParams/HotStartFrame in Bax.
    static const UInt Bax_ScanData_HotStartFrame = 0;
    // Default value of attribute /ScanData/AcqParams/LaserOnFrame in Bax.
    static const UInt Bax_ScanData_LaserOnFrame = 0;
    // Default value of attribute /ScanData/AcqParams/FrameRate in Bax.
    static const float Bax_ScanData_FrameRate = 80.047035;

    // Default value of attribute /ScanData/RunInfo/RunCode in Bax.
    static const std::string Bax_ScanData_RunCode = "Bam2Bax_Run_Code";
    // Default value of attribute /ScanData/DyeSet/BaseMap in Bax.
    static const std::string Bax_ScanData_BaseMap = PacBio::AttributeValues::ScanData::DyeSet::basemap;
    // Default value of attribute /Regions/RegionTypes in Bax.
    static const std::vector<std::string> Bax_Regions_RegionTypes = PacBio::AttributeValues::Regions::regiontypes;
}

template <class T_HDFWRITER>
class Bam2BaxConverter : public IConverter
{
public:
    Bam2BaxConverter(Settings & settings)
    :IConverter(settings) {}

    ~Bam2BaxConverter(void) {}

    bool Run(void) {return ConvertFile();}

protected:
    bool ConvertFile(void);
};

#include "Bam2BaxConverterImpl.hpp"
#endif
