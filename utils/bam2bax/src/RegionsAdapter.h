// Author: Yuan Li

#ifndef _REGIONS_ADAPTER_H_
#define _REGIONS_ADAPTER_H_

#include "pbdata/reads/RegionAnnotation.hpp"
#include "RegionTypeAdapter.h"

class RegionsAdapter {
public:
    /// \name \{
    /// Converts PacBio::BAM::VirtualRegion to RegionAnnotation
    /// in pbdata.
    /// VirtualRegion has four fields, including
    ///      * region type,
    ///      * region start,
    ///      * region end,
    ///      * region score.
    /// RegionAnnotation has five fields, including
    ///      * holeNumber,        --> missing in VirtualRegion
    ///      * region type index, --> not region type.
    ///      * region start,
    ///      * region end,
    ///      * region score
    /// \note region type has to be converted to region type index.
    /// \}
    /// \param[in] holeNumber, zmw hole number, which is absent in VirtualRegion,
    ///            but is required by RegionAnnotation
    /// \param[in] vr, virtual region, which contains region type, 
    ///            region start, region end and region score, out of which
    ///            region type has to be converted to region type index.
    /// \param[in] regionTypes, a table to look up region types according 
    ///            to region type index.
    /// \returns a RegionAnnotation object
    static inline
    RegionAnnotation ToRegionAnnotation(const UInt holeNumber, 
                                        const PacBio::BAM::VirtualRegion & vr, 
                                        std::vector<RegionType> & regionTypes);
    /*
    /// Create a RegionAnnotaion object.
    inline RegionAnnotation CreateRegionAnnotation(const UInt holeNumber, const PacBio::BAM::VirtualRegion & vr, std::vector<RegionType> & regionTypes);
    */

    /// Comparison between two RegionAnnotations to decide their order in H5 RegionTable.
    static inline 
    bool CmpRegionAnnotations(const RegionAnnotation & l, 
                              const RegionAnnotation & r);

    //inline std::vector<RegionAnnotation> RegionAnnotationsFromVirtualPolymeraseRead (const PacBio::BAM::VirtualPolymeraseBamRecord & record);
    
    /// Creates a vector of RegionAnnotations from a virtual polymerase bam record.
    /// \returns a vector of RegionAnnotations created from a virtual polymerase bam record.
    /// \param[in] record, input Virtual Polymerase Bam Record.
    /// \param[in] regionTypes, a table to look up region types according 
    ///            to region type index.
    static inline 
    std::vector<RegionAnnotation> ToRegionAnnotations(
            const PacBio::BAM::VirtualPolymeraseBamRecord & record, 
            std::vector<RegionType> & regionTypes);

    /// \}
}; // class RegionsAdapter


RegionAnnotation RegionsAdapter::ToRegionAnnotation(
    const UInt holeNumber, 
    const PacBio::BAM::VirtualRegion & vr, 
    std::vector<RegionType> & regionTypes) 
{
    int index = RegionTypeAdapter::ToRegionTypeIndex(vr.type, regionTypes);
    if (vr.type == PacBio::BAM::VirtualRegionType::HQREGION and 
        vr.beginPos == vr.endPos) {
        // bug 29935, by convention, use HQREGION 0, 0, 0 if no HQREGION is found.
        return RegionAnnotation(holeNumber, index, 0, 0, 0);
    } else return RegionAnnotation(holeNumber, index, vr.beginPos, vr.endPos, vr.score);
}

bool RegionsAdapter::CmpRegionAnnotations(const RegionAnnotation & l, 
                                          const RegionAnnotation & r) 
{
    assert(l.GetHoleNumber() == r.GetHoleNumber());
    if (l.GetTypeIndex() == r.GetTypeIndex()) {
        if (l.GetStart() == r.GetStart()) {
            return l.GetEnd() > r.GetEnd();
        } else {
            return l.GetStart() < r.GetStart();
        }
    } else {
        return (l.GetTypeIndex() < r.GetTypeIndex());
    }
}

std::vector<RegionAnnotation> RegionsAdapter::ToRegionAnnotations (
        const PacBio::BAM::VirtualPolymeraseBamRecord & record,
        std::vector<RegionType> & regionTypes) {

    auto virtualRegionMap_ = record.VirtualRegionsMap();
    std::vector<RegionAnnotation> ret;
    for (auto it = virtualRegionMap_.begin(); it != virtualRegionMap_.end(); it++) {
        for (PacBio::BAM::VirtualRegion vr: it->second) {
            if (RegionTypeAdapter::IsConvertible(vr.type, regionTypes)) {
                RegionAnnotation annotation = RegionsAdapter::ToRegionAnnotation(record.HoleNumber(), vr, regionTypes);
                if ((vr.type == PacBio::BAM::VirtualRegionType::ADAPTER or
                     vr.type == PacBio::BAM::VirtualRegionType::HQREGION) and 
                     record.HasReadAccuracy()) {
                    float rq = record.ReadAccuracy(); 
                    annotation.SetScore((rq >= 1.0)?(int(rq)):(int(rq * 1000)));
                }
                ret.push_back(annotation);
            } // Some region types such as LQRegion can not be converted
        }
    }
    std::sort(ret.begin(), ret.end(), RegionsAdapter::CmpRegionAnnotations); 

    return ret;
}

#endif
