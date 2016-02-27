// Author: Derek Barnett

#include "Settings.h"
#include "OptionParser.h"
#include <boost/algorithm/string.hpp>
#include <HDFNewBasReader.hpp>
#include <pbbam/DataSet.h>
#include <sstream>
using namespace std;

namespace internal {

static
vector<string> BaxFilenamesFromXml(const string& xmlFilename)
{
    using namespace PacBio::BAM;

    try {
        vector<string> filenames;

        DataSet dataset(xmlFilename);
        const vector<string> resources = dataset.ResolvedResourceIds();
        for (const string& resource : resources) {
            cerr << resource << endl;
            const boost::iterator_range<string::const_iterator> baxFound = boost::algorithm::ifind_first(resource, ".bax.h5");
            if (!baxFound.empty()) 
                filenames.push_back(resource);
        }
        return filenames;

    } catch (std::exception&) {
        // TODO: report error
        return vector<string>();
    }
}

static
vector<string> FilenamesFromFofn(const string& fileName)
{
    vector<string> retval;
    ifstream in_stream;
    string line;

    in_stream.open(fileName);

    while(!in_stream.eof())
    {
        in_stream >> line;
        if (!line.empty())
            retval.push_back(line);
        line.clear();
    }

    return retval;
}        

static
bool isBasH5(const string& fileName)
{
    return boost::ends_with(boost::to_lower_copy(fileName), ".bas.h5");
}

static
void H5FilenamesFromBasH5(const string& basFileName,
                          vector<string>* const output)
{
    HDFNewBasReader reader;
    if (reader.Initialize(basFileName))
        for (const auto& baxFileName : reader.GetBaxFileNames())
            output->push_back(baxFileName);
    else
        output->push_back(basFileName);
}

} // namespace internal

// option names
const char* Settings::Option::datasetXml_     = "datasetXml";
const char* Settings::Option::hqRegionMode_   = "hqRegionMode";
const char* Settings::Option::input_          = "input";
const char* Settings::Option::fofn_           = "fofn";
const char* Settings::Option::losslessFrames_ = "losslessFrames";
const char* Settings::Option::output_         = "output";
const char* Settings::Option::polymeraseMode_ = "polymeraseMode";
const char* Settings::Option::pulseFeatures_  = "pulseFeatures";
const char* Settings::Option::subreadMode_    = "subreadMode";
const char* Settings::Option::ccsMode_        = "ccsMode";
const char* Settings::Option::outputXml_      = "outputXml";
const char* Settings::Option::sequelPlatform_ = "sequelPlatform";

Settings::Settings(void)
    : mode(Settings::SubreadMode)
    , isSequelInput_(false)
    , usingDeletionQV(true)
    , usingDeletionTag(true)
    , usingInsertionQV(true)
    , usingIPD(true)
    , usingMergeQV(true)
    , usingPulseWidth(false)
    , usingSubstitutionQV(true)
    , usingSubstitutionTag(false)
    , losslessFrames(false)
{ }

Settings Settings::FromCommandLine(optparse::OptionParser& parser,
                                   int argc,
                                   char *argv[])
{
    Settings settings;

    // general program info
    settings.program = parser.prog();
    settings.description = parser.description();
    settings.version = parser.version();
    for (int i = 1; i < argc; ++i) {
        settings.args.append(argv[i]);
        settings.args.append(" ");
    }

    const optparse::Values options = parser.parse_args(argc, argv);

    // output prefix
    // TODO: output dir ??
    settings.outputBamPrefix = options[Settings::Option::output_];
    settings.outputXmlFilename = options[Settings::Option::outputXml_];

    // input files from dataset XML ?
    if ( options.is_set(Settings::Option::datasetXml_) ) {
        settings.datasetXmlFilename = options[Settings::Option::datasetXml_];
        settings.inputBaxFilenames = internal::BaxFilenamesFromXml(settings.datasetXmlFilename);
    }

    // input files from fofn ?
    else if ( options.is_set(Settings::Option::fofn_))
    {
        settings.fofnFilename = options[Settings::Option::fofn_];
        settings.inputFilenames = internal::FilenamesFromFofn(settings.fofnFilename);        
    }

    // else input files command-line args
    else
        settings.inputFilenames = parser.args();

    // Process input files to convert Bas.H5 --> Bax.h5 as needed
    for (const std::string& fn : settings.inputFilenames)
    {
        if (internal::isBasH5(fn))
            internal::H5FilenamesFromBasH5(fn, &settings.inputBaxFilenames);
        else
            settings.inputBaxFilenames.push_back(fn);
    }

    if (settings.inputBaxFilenames.empty())
        settings.errors_.push_back("missing input BAX files.");

    // mode
    const bool isSubreadMode =
            options.is_set(Settings::Option::subreadMode_) ? options.get(Settings::Option::subreadMode_)
                                                           : false;
    const bool isHQRegionMode =
            options.is_set(Settings::Option::hqRegionMode_) ? options.get(Settings::Option::hqRegionMode_)
                                                            : false;
    const bool isPolymeraseMode =
            options.is_set(Settings::Option::polymeraseMode_) ? options.get(Settings::Option::polymeraseMode_)
                                                              : false;
    const bool isCCS =
            options.is_set(Settings::Option::ccsMode_) ? options.get(Settings::Option::ccsMode_)
                                                       : false;

    int modeCount = 0;
    if (isSubreadMode)    ++modeCount;
    if (isHQRegionMode)   ++modeCount;
    if (isPolymeraseMode) ++modeCount;
    if (isCCS)            ++modeCount;

    if (modeCount == 0)
        settings.mode = Settings::SubreadMode;
    else if (modeCount == 1) {
        if (isSubreadMode)    settings.mode = Settings::SubreadMode;
        if (isHQRegionMode)   settings.mode = Settings::HQRegionMode;
        if (isPolymeraseMode) settings.mode = Settings::PolymeraseMode;
        if (isCCS)            settings.mode = Settings::CCSMode;
    }
    else
        settings.errors_.push_back("multiple modes selected");

    // platform
    settings.isSequelInput_ = options.is_set(Settings::Option::sequelPlatform_) ? options.get(Settings::Option::sequelPlatform_)
                                                                                : false;

    // frame data encoding
    settings.losslessFrames = options.is_set(Settings::Option::losslessFrames_) ? options.get(Settings::Option::losslessFrames_)
                                                                                : false;

    // pulse features list
    if (options.is_set(Settings::Option::pulseFeatures_)) {

        // ignore defaults
        settings.usingDeletionQV = false;
        settings.usingDeletionTag = false;
        settings.usingInsertionQV = false;
        settings.usingIPD = false;
        settings.usingMergeQV = false;
        settings.usingPulseWidth = false;
        settings.usingSubstitutionQV = false;
        settings.usingSubstitutionTag = false;

        // apply user-requested features
        stringstream stream(options[Settings::Option::pulseFeatures_]);
        string feature;
        while(std::getline(stream, feature, ',')) {
            if      (feature == "DeletionQV")      settings.usingDeletionQV = true;
            else if (feature == "DeletionTag")     settings.usingDeletionTag = true;
            else if (feature == "InsertionQV")     settings.usingInsertionQV = true;
            else if (feature == "IPD")             settings.usingIPD = true;
            else if (feature == "MergeQV")         settings.usingMergeQV = true;
            else if (feature == "PulseWidth")      settings.usingPulseWidth = true;
            else if (feature == "SubstitutionQV")  settings.usingSubstitutionQV = true;
            else if (feature == "SubstitutionTag") settings.usingSubstitutionTag = true;
            else
                settings.errors_.push_back(string("unknown pulse feature: ") + feature);
        }
    }

#ifdef DEBUG_SETTINGS

    string modeString;
    if (settings.mode == Settings::SubreadMode)
        modeString = "subread";
    else if (settings.mode == Settings::HQRegionMode)
        modeString = "hqRegion";
    else if (settings.mode == Settings::PolymeraseMode)
        modeString = "polymerase";
    else
        modeString = "ccs";

    string platformString = settings.isSequelInput_ ? "Sequel" : "RS";

    cerr << "CommandLine: " << settings.program << " " << settings.args << endl
         << "Description: " << settings.description << endl
         << "Version:     " << settings.version << endl
         << "Mode:        " << modeString << endl
         << "Platform:    " << platformString << endl
         << "DeletionQV?:      " << ( settings.usingDeletionQV ? "yes" : "no" ) << endl
         << "DeletionTag?:     " << ( settings.usingDeletionTag ? "yes" : "no" ) << endl
         << "InsertionQV?:     " << ( settings.usingInsertionQV ? "yes" : "no" ) << endl
         << "IPD?:             " << ( settings.usingMergeQV ? "yes" : "no" ) << endl
         << "MergeQV?:         " << ( settings.usingIPD ? "yes" : "no" ) << endl
         << "PulseWidth?:      " << ( settings.usingPulseWidth ? "yes" : "no" ) << endl
         << "SubstitutionQV?:  " << ( settings.usingSubstitutionQV ? "yes" : "no" ) << endl
         << "SubstitutionTag?: " << ( settings.usingSubstitutionTag ? "yes" : "no" ) << endl;
#endif

    return settings;
}
