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


    // Write metadata.xml to parent directory of Bax.h5.
    if (not settings_.outputMetadataFilename.empty())
        MetadataWriter metaWriter_(settings_.outputMetadataFilename, 
                                   rg, 
                                   settings_.outputAnalysisDirname);

    T_HDFWRITER writer(outfn, 
            rg.BasecallerVersion(), 
            scandata.BaseMap(),
            qvs,
            Bam2BaxDefaults::Bax_Regions_RegionTypes);

    if (settings_.traceFilename.empty()) {
        writer.WriteScanData(scandata);
    } else {
        HDFFile traceFile;
        traceFile.Open(settings_.traceFilename, H5F_ACC_RDONLY);
        writer.CopyObject(traceFile, "/ScanData"); 
        traceFile.Close();
    }
       
    if (not settings_.subreadsBamFilename.empty() and 
        not settings_.scrapsBamFilename.empty()) {

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
    }
            
    return errors_.empty();
}
#endif
