// Author: Yuan Li
#ifndef BAM2BAX_ICONVERTER_H_
#define BAM2BAX_ICONVERTER_H_

#include <string>
#include <vector>
#include <algorithm>
#include "pbdata/Enumerations.h"
#include "pbbam/BamFile.h"
#include "pbbam/BamHeader.h"
#include "pbbam/ReadGroupInfo.h"
#include "pbbam/virtual/VirtualPolymeraseReader.h"
#include "pbbam/virtual/VirtualPolymeraseBamRecord.h"
#include "pbbam/virtual/VirtualRegion.h"
#include "pbbam/virtual/VirtualRegionType.h"
#include "pbbam/virtual/VirtualRegionTypeMap.h"
#include "HDFWriterBase.hpp"
#include "HDFBaxWriter.hpp"
#include "HDFPulseWriter.hpp"
#include "RegionsAdapter.h"
#include "Settings.h"
#include "Bam2BaxInternal.h"

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

class Converter {
public:
    Converter(Settings & settings);
    ~Converter(void);

public:
    std::vector<std::string> Errors(void) const;
    bool Run();

protected:
    void AddErrorMessage(const std::string & errmsg) {
        errors_.push_back(errmsg);
    }

protected:
    // protected variables
    Settings& settings_;
    ScanData* scanData_;
    HDFWriterBase* writer_;
    PacBio::BAM::BamFile* bamfile_;
    std::vector<std::string> errors_;

private:
    void MockScanData(PacBio::BAM::ReadGroupInfo& rg);
    void InitializeWriter(const std::string& bcvers, 
                          const std::vector<PacBio::BAM::BaseFeature>& qvs);
};
#endif
