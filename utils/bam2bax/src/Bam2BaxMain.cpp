// Author: Yuan Li

#include "Bam2Bax.h"
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
    parser.description("bam2bax converts the PacBio BAM format into bax.h5 format.");
    parser.prog("bam2bax");
    parser.version("1.0.0.170337");
    parser.add_version_option(true);
    parser.add_help_option(true);

    auto ioGroup = optparse::OptionGroup(parser, "Input/output files");
    ioGroup.add_option("")
           .dest(Settings::Option::input_)
	       .metavar("movie.subreads.bam movie.scraps.bam | movie.polymerase.bam") 
           .help("Input a movie.polymerase.bam. Or a movie.subreads.bam and a movie.scraps.bam");
    ioGroup.add_option("-o")
           .dest(Settings::Option::output_)
	       .metavar("STRING")
           .help("Prefix of output filenames. Movie name will be used if no prefix provided");
    ioGroup.add_option("--metadata")
           .dest(Settings::Option::metadata_)
           .action("store_true")
           .help("Write metadata.xml to the upper directory of output file.");
    parser.add_option_group(ioGroup);

    auto modeGroup = optparse::OptionGroup(parser, "Output file types (mutually exclusive:)");
    modeGroup.add_option("--base")
             .dest(Settings::Option::baseMode_)
             .metavar("")
             .action("store_true")
             .help("Output bax.h5 (default)");
    modeGroup.add_option("--pulse")
             .dest(Settings::Option::pulseMode_)
             .metavar("")
             .action("store_true")
             .help("Output pls.h5");
    modeGroup.add_option("--baseMap")
             .dest(Settings::Option::baseMap_)
             .metavar(Settings::OptionValue::baseMap_)
             .help("Set /ScanData/DyeSet/BaseMap, mapping channels to bases.");
    modeGroup.add_option("--ignoreQV")
             .dest(Settings::Option::ignoreQV_)
             .metavar("")
             .action("store_true")
             .help("Don't save QVs in ouptut file.");
    parser.add_option_group(modeGroup);

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
    return Bam2Bax::Run(settings);
}
