#ifndef _TEST_CONSTANTS_
#define _TEST_CONSTANTS_
#include "Enumerations.h"
#include <string>
#include <vector>
#include <map>
#include "SMRTSequence.hpp"
#include "reads/ScanData.hpp"

#include "HDFFile.hpp"
#include "HDFScanDataWriter.hpp"
#include "HDFBaseCallsWriter.hpp"
#include "HDFBaxWriter.hpp"
#include "HDFBasReader.hpp"
#include "HDFScanDataReader.hpp"

namespace tests {
    // Setup variables for test HDFScanDataWriter.
    const PlatformId platformid = Springfield;
    const std::string movieName = "thisismymovie";
    const std::string runCode = "thisismyruncode";
    const std::string whenStarted = "2015-01-05";
    const float frameRate = 75.0;
    const unsigned int numFrames = 539900;
    const std::map<char, size_t> baseMap = {{'T', 0}, {'G', 1}, {'A', 2}, {'C', 3}};
    const std::map<char, size_t> randomBaseMap = {{'A', 1}, {'C', 3}, {'G', 2}, {'T', 0}};
    const std::string bindingKit = "thisismybindingkit";
    const std::string sequencingKit = "thisismysequencingkit";

    // Setup variables for test HDFBaseCallersWriter.
    const std::vector<std::string> QVNames = 
        {"DeletionTag", "DeletionQV", 
         "InsertionQV", "MergeQV", 
         "SubstitutionQV", "SubstitutionTag", 
         "PreBaseFrames", "WidthInFrames",
         "HQRegionSNR", "ReadScore"};

    const std::vector<PacBio::BAM::BaseFeature> QVEnums =
        {   PacBio::BAM::BaseFeature::DELETION_TAG
          , PacBio::BAM::BaseFeature::DELETION_QV
          , PacBio::BAM::BaseFeature::INSERTION_QV
          , PacBio::BAM::BaseFeature::MERGE_QV
          , PacBio::BAM::BaseFeature::SUBSTITUTION_QV
          , PacBio::BAM::BaseFeature::SUBSTITUTION_TAG
          , PacBio::BAM::BaseFeature::IPD
          , PacBio::BAM::BaseFeature::PULSE_WIDTH };


    const unsigned int len = 2;
    const unsigned int holeNumber = 110;
    const unsigned char holeStatus = '8';
    const int readScoreInt = 760;
    const float snra = 0.11;
    const float snrc = 0.22;
    const float snrg = 0.33;
    const float snrt = 0.44;

    const std::string myseq = "AT";
    const std::string deletionTag = "NG";
    const std::string deletionQV = "&2";
    const std::string insertionQV = "+%";
    const std::string mergeQV = "(.";
    const std::string substitutionQV = "/<";
    const std::string substitutionTag = "NC";
    const std::string preBaseFrames = "C0";
    const std::string widthInFrames = "8;";

    const std::string PULSEDATA = "PulseData";

    const std::string basecallerVersion = "2.0.1.354";

    template <typename T>
    inline void SetData(T dst, const std::string & src, const unsigned int len) {
        for (unsigned int i = 0; i < len; i++) 
            dst[i] = src.c_str()[i];
    }

    template <typename T>
    inline bool CmpData(T l, T r, const unsigned int len) {
        for (unsigned int i = 0; i < len; i++) {
            if (l[i] != r[i]) return false;
        }
        return true;
    }

    inline void make_scandata(ScanData & scandata, 
                              const std::map<char, size_t> & baseMap_) {
    scandata.PlatformID(platformid)
        .MovieName(tests::movieName)
        .WhenStarted(tests::whenStarted)
        .RunCode(runCode)
        .NumFrames(tests::numFrames)
        .FrameRate(tests::frameRate)
        .SequencingKit(tests::sequencingKit)
        .BindingKit(tests::bindingKit)
        .BaseMap(baseMap_);
    }

    inline void  make_smrtseq(SMRTSequence & seq) {
        // Create a SMRTSequence
        seq.Allocate(tests::len);

        memcpy(seq.seq, myseq.c_str(), tests::len * sizeof(char));

        SetData<QualityValueVector<unsigned char>>(seq.deletionQV, tests::deletionQV,    tests::len);
        SetData<Nucleotide *>(seq.deletionTag, tests::deletionTag,   tests::len);
        SetData<QualityValueVector<unsigned char>>(seq.insertionQV, tests::insertionQV,   tests::len);
        SetData<QualityValueVector<unsigned char>>(seq.mergeQV,     tests::mergeQV,       tests::len);
        SetData<QualityValueVector<unsigned char>>(seq.substitutionQV,  tests::substitutionQV,  tests::len);
        SetData<Nucleotide *>(seq.substitutionTag, tests::substitutionTag, tests::len);
        SetData<HalfWord *>(seq.preBaseFrames,    tests::preBaseFrames, tests::len);
        SetData<HalfWord *>(seq.widthInFrames,    tests::widthInFrames, tests::len);

        seq.zmwData.holeNumber = tests::holeNumber;
        seq.zmwData.holeStatus = tests::holeStatus;

        seq.readScore = tests::readScoreInt / 1000.0;

        seq.HQRegionSnr('A', tests::snra);
        seq.HQRegionSnr('C', tests::snrc);
        seq.HQRegionSnr('G', tests::snrg);
        seq.HQRegionSnr('T', tests::snrt);
    }

    // Write seq to outfn
    inline bool write_to(const std::string outfn, const ScanData & sd, const SMRTSequence & seq) {
        std::unique_ptr<HDFBaxWriter> writer;
        writer.reset(new HDFBaxWriter(outfn, tests::basecallerVersion, sd.BaseMap(), tests::QVEnums));
        bool ret = writer->WriteOneZmw(seq);
        writer.reset();
        return ret;
    }

    // Write seq to outfn
    inline bool write_to(const std::string outfn, const SMRTSequence & seq) {
        HDFFile outfile;
        outfile.Open(outfn, H5F_ACC_TRUNC);

        HDFGroup pulseDataGroup;
        outfile.rootGroup.AddGroup(tests::PULSEDATA);
        pulseDataGroup.Initialize(outfile.rootGroup, tests::PULSEDATA);

        std::unique_ptr<HDFBaseCallsWriter> writer;
        writer.reset(new HDFBaseCallsWriter(outfn, pulseDataGroup, tests::baseMap, tests::basecallerVersion, tests::QVEnums)); 

        bool ret = writer->WriteOneZmw(seq);
        writer.reset();
        outfile.Close();
        return ret;
    }

    // Read seq from fn
    inline int read_from(const std::string fn, SMRTSequence & seq) {
        // read the sequence from fn
        HDFBasReader reader;
        std::vector<std::string> qvn = tests::QVNames;
        reader.InitializeFields(qvn);
        reader.Initialize(fn);
        int count = 0;
        while(reader.GetNext(seq)) {
            count += 1;
        }
        return count;
    }

    inline void write_to(const std::string & outfn, ScanData & scandata) {
        HDFFile outfile;
        outfile.Open(outfn, H5F_ACC_TRUNC);
        HDFScanDataWriter writer(outfile.rootGroup);
        writer.Write(scandata);
        writer.Close();
        outfile.Close();
    }

    inline bool read_from(const std::string & infn, ScanData & scandata) {
        //read
        HDFScanDataReader reader;
        HDFFile infile;
        infile.Open(infn, H5F_ACC_RDONLY);
        HDFGroup * rootGroupPtr = &infile.rootGroup;

        if (not rootGroupPtr->ContainsObject("ScanData")) return false;
        
        reader.Initialize(rootGroupPtr);

        bool OK = reader.Read(scandata);
        reader.Close();
        infile.Close();
        return OK;
    }

};

#endif
