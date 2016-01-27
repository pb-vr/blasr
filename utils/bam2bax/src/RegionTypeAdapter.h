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

#ifndef _REGIONTYPE_ADAPTER_H_
#define _REGIONTYPE_ADAPTER_H_

#include <string>
#include <vector>

class RegionTypeAdapter {
public:
    /// \name \{
    /// Converts between RegionType and VirtualRegionType
    /// \returns true if input PacBio::BAM::VirtualRegionType object vrt can be converted to any RegionType in regionTypes.
    /// \param[in]  vrt, VirtualRegionType to be converted.
    /// \param[in]  regionTypes, valid RegionTypes which can be converted to.
    static inline 
    //bool IsConvertibleVirtualRegionType(PacBio::BAM::VirtualRegionType vrt, std::vector<RegionType> & regionTypes); 
    bool IsConvertible(PacBio::BAM::VirtualRegionType vrt, std::vector<RegionType> & regionTypes); 

    /// Converts PacBio::BAM::VirtualRegionType vrt to enum RegionType defined in pbdata/Enumeration.h
    //inline RegionType VirtualRegionTypeToRegionType(PacBio::BAM::VirtualRegionType vrt);
    static inline
    RegionType ToRegionType(PacBio::BAM::VirtualRegionType vrt);

    static inline
    RegionType ToRegionType(const std::string & type);

    static inline
    std::vector<RegionType> ToRegionTypes(const std::vector<std::string> & typeStrs);

    static inline
    RegionType ToRegionTypes(const std::string & str);

    static inline
    PacBio::BAM::VirtualRegionType ToVirtualRegionType(const std::string & str);

    /// Converts RegionType to PacBio::BAM::VirtualRegionType
    //inline PacBio::BAM::VirtualRegionType RegionTypeToVirtualRegionType(RegionType rt);
    static inline
    PacBio::BAM::VirtualRegionType ToVirtualRegionType(RegionType rt);
    
    /// Converts VirtualRegionType to RegionTypeIndex in regionTypes.
    /// \returns index of this region type in given regionTypes.
    static inline 
    int ToRegionTypeIndex(PacBio::BAM::VirtualRegionType vrt, std::vector<RegionType> & regionTypes);

    /// \}
}; // RegionTypeAdapter


/// Convert a BAM::VirtualRegionType to pbdata RegionAnnotataion.TypeIndex.
int RegionTypeAdapter::ToRegionTypeIndex(PacBio::BAM::VirtualRegionType vrt,
                                         std::vector<RegionType> & regionTypes) {
    RegionType rt =  ToRegionType(vrt); 
    std::vector<RegionType>::iterator it =  std::find(regionTypes.begin(), 
                                                      regionTypes.end(), 
                                                      rt);
    size_t index =  std::distance(regionTypes.begin(), it);
    assert(index !=  regionTypes.size());
    return static_cast<int>(index);
}

bool RegionTypeAdapter::IsConvertible(PacBio::BAM::VirtualRegionType vrt, 
                                      std::vector<RegionType> & regionTypes) 
{
    RegionType rt = ToRegionType(vrt);
    if (rt == UnknownRegionType) return false;
    std::vector<RegionType>::iterator it = std::find(regionTypes.begin(), regionTypes.end(), rt);
    return (it != regionTypes.end());
}

RegionType RegionTypeAdapter::ToRegionType(PacBio::BAM::VirtualRegionType vrt)
{ 
    if (vrt == PacBio::BAM::VirtualRegionType::SUBREAD)
        return Insert;
    else if (vrt == PacBio::BAM::VirtualRegionType::ADAPTER)
        return Adapter;
    else if (vrt == PacBio::BAM::VirtualRegionType::HQREGION)
        return HQRegion;
    else if (vrt == PacBio::BAM::VirtualRegionType::BARCODE)
        return BarCode;
    else 
        return UnknownRegionType;
    //e.g., No LQRegion defined in pbdata/Enumeration.h
}

RegionType RegionTypeAdapter::ToRegionType(const std::string & str) {
    std::string u_str = str;
    std::transform(u_str.begin(), u_str.end(), u_str.begin(), ::toupper);
    if (u_str == "INSERT" || u_str == "SUBREAD") {
        return Insert;
    } else if (u_str == "ADAPTER") {
        return Adapter;
    } else if (u_str == "HQREGION") {
        return HQRegion;
    } else if (u_str == "BARCODE") {
        return BarCode;
    } else {
        return UnknownRegionType;
    }
}

std::vector<RegionType> RegionTypeAdapter::ToRegionTypes(const std::vector<std::string> & typeStrs) {
    std::vector<RegionType> ret;
    for(auto str: typeStrs) 
        ret.push_back(ToRegionType(str));
    return ret;
}

PacBio::BAM::VirtualRegionType RegionTypeAdapter::ToVirtualRegionType(const std::string & str) {
    return ToVirtualRegionType(ToRegionType(str));
}

PacBio::BAM::VirtualRegionType RegionTypeAdapter::ToVirtualRegionType(RegionType rt) {
    if (rt == Insert) 
        return PacBio::BAM::VirtualRegionType::SUBREAD;
    else if (rt == Adapter)
        return PacBio::BAM::VirtualRegionType::ADAPTER;
    else if (rt == HQRegion)
        return PacBio::BAM::VirtualRegionType::HQREGION;
    else if (rt == BarCode)
        return PacBio::BAM::VirtualRegionType::BARCODE;
    else
        assert("Unable to convert RegionType to VirtualRegionType." == NULL);
}



#endif
