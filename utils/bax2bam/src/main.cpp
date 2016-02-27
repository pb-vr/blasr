// Author: Derek Barnett

#include "Bax2Bam.h"
#include "OptionParser.h"
#include "Settings.h"
#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

int main(int argc, char* argv[])
{
    // setup help & options
    optparse::OptionParser parser;
    parser.description("bax2bam converts the legacy PacBio basecall format (bax.h5) into the BAM basecall format.");
    parser.prog("bax2bam");
    parser.version("0.0.3");
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto ioGroup = optparse::OptionGroup(parser, "Input/output files");
    ioGroup.add_option("")
           .dest(Settings::Option::input_)
	   .metavar("movie.1.bax.h5 movie.2.bax.h5 ...")
           .help("Input files which should be from the same movie");
    ioGroup.add_option("--xml")
           .dest(Settings::Option::datasetXml_)
           .metavar("STRING")
           .help("DataSet XML file containing a list of movie names");
    ioGroup.add_option("-f", "--fofn")
           .dest(Settings::Option::fofn_)
           .metavar("STRING")
           .help("File-of-file-names containing a list of input files");
    ioGroup.add_option("-o")
           .dest(Settings::Option::output_)
	   .metavar("STRING")
           .help("Prefix of output filenames. Movie name will be used if no prefix provided");
    ioGroup.add_option("--output-xml")
           .dest(Settings::Option::outputXml_)
           .metavar("STRING")
           .help("Explicit output XML name. If none provided via this arg, bax2bam will use -o prefix (<prefix>.dataset.xml). "
                 "If that is not specified either, the output XML filename will be <moviename>.dataset.xml");
    parser.add_option_group(ioGroup);

    auto platformGroup = optparse::OptionGroup(parser, "Input sequencing platform");
    platformGroup.add_option("--sequel-input")
                 .dest(Settings::Option::sequelPlatform_)
                 .action("store_true")
                 .help("Specify that input data is from Sequel. "
                       "bax2bam will assume RS unless this option is specified");

    auto modeGroup = optparse::OptionGroup(parser, "Output read types (mutually exclusive:)");
    modeGroup.add_option("--subread")
             .dest(Settings::Option::subreadMode_)
             .action("store_true")
             .help("Output subreads (default)");
    modeGroup.add_option("--hqregion")
             .dest(Settings::Option::hqRegionMode_)
             .action("store_true")
             .help("Output HQ regions");
    modeGroup.add_option("--polymeraseread")
             .dest(Settings::Option::polymeraseMode_)
             .action("store_true")
             .help("Output full polymerase read");
    modeGroup.add_option("--ccs")
             .dest(Settings::Option::ccsMode_)
             .action("store_true")
             .help("Output CCS sequences");
    parser.add_option_group(modeGroup);

    auto featureGroup = optparse::OptionGroup(parser, "Pulse feature options");
    featureGroup.group_description("Configure pulse features in the output BAM. Supported features include:\n"
                                   "    Pulse Feature:    BAM tag:  Default:\n"
                                   "    DeletionQV        dq        Y\n"
                                   "    DeletionTag       dt        Y\n"
                                   "    InsertionQV       iq        Y\n"
                                   "    IPD               ip        Y\n"
                                   "    PulseWidth        pw        N\n"
                                   "    MergeQV           mq        Y\n"
                                   "    SubstitutionQV    sq        Y\n"
                                   "    SubstitutionTag   st        N\n"
                                   "If this option is used, then only those features listed will be included, "
                                   "regardless of the default state."
                                   );
    featureGroup.add_option("--pulsefeatures")
                .dest(Settings::Option::pulseFeatures_)
                .metavar("STRING")
                .help("Comma-separated list of desired pulse features, using the names listed above.\n");
    featureGroup.add_option("--losslessframes")
                .dest(Settings::Option::losslessFrames_)
                .action("store_true")
                .help("Store full, 16-bit IPD/PulseWidth data, instead of (default) downsampled, 8-bit encoding.");
    parser.add_option_group(featureGroup);

    // parse command line
    Settings settings = Settings::FromCommandLine(parser, argc, argv);
    if (!settings.errors_.empty()) {
        cerr << endl;
        for (const auto e : settings.errors_)
            cerr << "ERROR: " << e << endl;
        cerr << endl;
        parser.print_help();
        return EXIT_FAILURE;
    }

    // main conversion
    return Bax2Bam::Run(settings);
}
