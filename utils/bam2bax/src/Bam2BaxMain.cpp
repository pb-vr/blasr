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
