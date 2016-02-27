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
