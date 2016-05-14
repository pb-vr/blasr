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
#include "HDFFile.hpp"
#include "RegionsAdapter.h"
#include "IConverter.h"



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
