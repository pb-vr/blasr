#include "IConverter.h"

IConverter::IConverter(Settings & settings)
:settings_(settings) { 
    writer_ = NULL;
    scanData_ = NULL;

    std::string infn = settings_.subreadsBamFilename;

    if (infn.empty()) infn = settings_.polymeraseBamFilename;

    PacBio::BAM::BamFile bamfile(infn);
    PacBio::BAM::BamHeader bamheader = bamfile.Header();

    if (bamheader.ReadGroups().size() != 1) {
        AddErrorMessage("Bam file must contain reads from exactly one SMRTCell.");
        // XXX: Throw initialization exception
    }
    PacBio::BAM::ReadGroupInfo rg = bamheader.ReadGroups()[0];
    MockScanData(rg);

    // FIXME: pbbam needs to provide an API which returns BaseFeatures in read group
    std::vector<PacBio::BAM::BaseFeature> qvs = settings_.ignoreQV ? std::vector<PacBio::BAM::BaseFeature>({}) : internal::QVEnumsInFirstRecord(bamfile);

    // Regions attribute RegionTypes, which defines supported region types in ORDER.
    std::vector<RegionType> regionTypes = RegionTypeAdapter::ToRegionTypes(Bam2BaxDefaults::Bax_Regions_RegionTypes);

    InitializeWriter(rg.BasecallerVersion(), qvs);
}

IConverter::~IConverter(void) {
    if (scanData_ != NULL) delete scanData_;
    if (writer_ != NULL) delete writer_;
}

std::vector<std::string> IConverter::Errors(void) const {
    return errors_;
}

bool IConverter::Run() {
    if (settings_.traceFilename.empty()) {
        writer_->WriteScanData(*scanData_);
    } else {
        HDFFile traceFile;
        traceFile.Open(settings_.traceFilename, H5F_ACC_RDONLY);
        writer_->CopyObject(traceFile, "/ScanData"); 
        // XXX: setup the inverse gain if writing pulses
        traceFile.Close();
    }

    return errors_.empty();
}

void IConverter::MockScanData(PacBio::BAM::ReadGroupInfo& rg) {
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

void IConverter::InitializeWriter(const std::string& bcvers, 
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
        // XXX: Throw initialization exception
    }
}
