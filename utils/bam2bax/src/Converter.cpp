#include "Converter.h"

Converter::Converter(Settings const& settings)
:settings_(settings) { 
    writer_ = NULL;
    scanData_ = NULL;

    std::string infn = settings_.subreadsBamFilename;

    if (infn.empty()) infn = settings_.polymeraseBamFilename;

    bamfile_ = new PacBio::BAM::BamFile(infn);
    PacBio::BAM::BamHeader bamheader = bamfile_->Header();

    if (bamheader.ReadGroups().size() != 1) {
        AddErrorMessage("Bam file must contain reads from exactly one SMRTCell.");
        // XXX: Throw initialization exception
    }
    PacBio::BAM::ReadGroupInfo rg = bamheader.ReadGroups()[0];
    MockScanData(rg);

    // FIXME: pbbam needs to provide an API which returns BaseFeatures in read group
    std::vector<PacBio::BAM::BaseFeature> qvs = settings_.ignoreQV ? std::vector<PacBio::BAM::BaseFeature>({}) : internal::QVEnumsInFirstRecord(*bamfile_);

    InitializeWriter(rg.BasecallerVersion(), qvs);
}

Converter::~Converter(void) {
    if (scanData_ != NULL) delete scanData_;
    if (writer_ != NULL) delete writer_;
    delete bamfile_;
}

std::vector<std::string> Converter::Errors(void) const {
    return errors_;
}

bool Converter::Run() {
    if (settings_.traceFilename.empty()) {
        writer_->WriteScanData(*scanData_);
    } else {
        HDFFile traceFile;
        traceFile.Open(settings_.traceFilename, H5F_ACC_RDONLY);
        writer_->CopyObject(traceFile, "/ScanData"); 
        if (settings_.mode == Settings::PulseMode) {
            SetInverseGain(traceFile);
        }
        traceFile.Close();
    }

    // Regions attribute RegionTypes, which defines supported region types in ORDER.
    std::vector<RegionType> regionTypes = 
        RegionTypeAdapter::ToRegionTypes(Bam2BaxDefaults::Bax_Regions_RegionTypes);

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
            if (not writer_->WriteOneZmw(smrt, ras) or not writer_->Errors().empty()) { break; }
            writer_->Flush();
        }
        if (not settings_.ignoreQV) writer_->WriteFakeDataSets();
        for (auto error: writer_->Errors()) { AddErrorMessage(error); }
    } else if (not settings_.polymeraseBamFilename.empty()) {
        // Read polymerase reads from polymerase.bam directly.
        PacBio::BAM::EntireFileQuery query(*bamfile_);
        for (auto record: query) {
            SMRTSequence smrt;
            smrt.Copy(record, true);
            RegionAnnotation ra(record.HoleNumber(), 
                                RegionTypeAdapter::ToRegionTypeIndex(PacBio::BAM::VirtualRegionType::HQREGION, regionTypes),
                                0, 0, 0);
            std::vector<RegionAnnotation> ras({ra});
            if (not writer_->WriteOneZmw(smrt, ras) or not writer_->Errors().empty()) { break; }
        }
        if (not settings_.ignoreQV) writer_->WriteFakeDataSets();
        for (auto error: writer_->Errors()) { AddErrorMessage(error); }
    }

    return errors_.empty();
}

void Converter::MockScanData(PacBio::BAM::ReadGroupInfo& rg) {
    // Construct AcqParams
    AcqParams acqParams(Bam2BaxDefaults::Bax_ScanData_AduGain,
                        Bam2BaxDefaults::Bax_ScanData_CameraGain,
                        Bam2BaxDefaults::Bax_ScanData_CameraType,
                        Bam2BaxDefaults::Bax_ScanData_HotStartFrame,
                        Bam2BaxDefaults::Bax_ScanData_LaserOnFrame);

    // Construct scandata.
    scanData_ = new ScanData(acqParams);
    scanData_->PlatformID(Sequel) // assume sequel movie 
             .MovieName(rg.MovieName()) // should be reliable now
             .WhenStarted(rg.Date())
             .RunCode(Bam2BaxDefaults::Bax_ScanData_RunCode)  // bam does not contain RunCode
             .NumFrames(Bam2BaxDefaults::Bax_ScanData_NumFrames) // bam does not contain NumFrames
             .FrameRate(Bam2BaxDefaults::Bax_ScanData_FrameRate) // Ignore bam header FrameRate.
             .SequencingKit(rg.SequencingKit())
             .BindingKit(rg.BindingKit())
             .BaseMap(settings_.baseMap);
}

void Converter::InitializeWriter(const std::string& bcvers, 
        const std::vector<PacBio::BAM::BaseFeature>& qvs) 
    {
    std::string outfn = settings_.outputBaxFilename;
    Settings::Mode mode = settings_.mode;
    
    if (mode == Settings::BaseMode) {
        std::cout << "Converting BAM to bax.h5." << std::endl;
        writer_ = new HDFBaxWriter(outfn, bcvers,
            scanData_->BaseMap(), qvs, Bam2BaxDefaults::Bax_Regions_RegionTypes);
    } else if (mode == Settings::PulseMode) {
        std::cout << "Converting BAM to plx.h5." << std::endl;
        writer_ = new HDFPulseWriter(outfn, bcvers,
            scanData_->BaseMap(), qvs, Bam2BaxDefaults::Bax_Regions_RegionTypes);
    } else {
        std::cerr << "UNKNOWN mode." << settings_.mode << std::endl;
        throw std::exception();
    }
}

void Converter::SetInverseGain(HDFFile& traceFile) {
    H5::Group acqGrp = traceFile.hdfFile.openGroup("/ScanData/AcqParams");
    H5::Attribute aduAttr = acqGrp.openAttribute("AduGain");
    float igain;
    H5::DataType* dt = new H5::DataType(H5::PredType::IEEE_F32LE);
    aduAttr.read(*dt, &igain);
    HDFPulseWriter* pw = static_cast<HDFPulseWriter*>(writer_);
    pw->SetInverseGain(igain);
}
