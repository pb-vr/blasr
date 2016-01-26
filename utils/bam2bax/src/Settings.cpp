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

#include "Settings.h"
#include "OptionParser.h"
#include "StringUtils.hpp"
#include <unistd.h> // getcwd
#include <string>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include "PacBioDefs.h"

#define DEBUG_SETTINGS

using namespace std;

namespace internal {

static inline
bool StartsWith(const string& input, const string& query)
{ return input.find(query) != string::npos; }

static
std::string GetMovienameFromFilename(const std::string & filename) {
    std::vector<std::string> tokens; 
    Splice(filename, "/", tokens);
    std::string tmp = tokens.back();
    Splice(tmp, ".", tokens);
    return tokens.front();
}

static 
bool IsAbsolutePath(const std::string & file) 
{
    return (file.find("/") == 0);
}

static
std::string CurrentWorkingDirectory(void) 
{
    char result[FILENAME_MAX] = { };
    if (getcwd(result, FILENAME_MAX) == nullptr)
        return std::string();
    return std::string(result);
}

static 
std::string DirectoryPath(const std::string & file) 
{
    // Return either relative or absolute directory path of file
    std::size_t pos = file.rfind('/');
    if (IsAbsolutePath(file)) {
        if (pos != std::string::npos) 
            return file.substr(0, pos);
        else 
            return std::string();
    } else {
        if (pos != std::string::npos)
            return CurrentWorkingDirectory() + "/" + file.substr(0, pos);
        else
            return CurrentWorkingDirectory();
    }
}

static 
std::string EXEC(const std::string & cmd) {
    // execute a cmd and return output
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}

static  // return parent directory name, no absolute path
std::string ParentDirectoryName (const std::string & folder) {
    std::stringstream ss;
    ss << "realpath " << folder << " | xargs basename";
    return EXEC(ss.str());
}

} // namespace internal

const char* Settings::Option::input_        = "input";
const char* Settings::Option::output_       = "output";
const char* Settings::Option::metadata_     = "metadata";
const char* Settings::Option::baseMode_     = "base";
const char* Settings::Option::pulseMode_    = "pulse";
const char* Settings::Option::baseMap_      = "basemap";
const char* Settings::Option::ignoreQV_     = "ignoreQV";
const char* Settings::OptionValue::baseMap_ = PacBio::AttributeValues::ScanData::DyeSet::basemap.c_str();

Settings::Settings(void)
    : mode(Settings::BaseMode)
    , ignoreQV(false)
    , baseMap(Settings::OptionValue::baseMap_)
{}

Settings Settings::FromCommandLine(optparse::OptionParser& parser,
                                   int argc,
                                   char *argv[],
                                   bool forcePulseMode)
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

    // mode
    settings.ignoreQV = options.is_set(Settings::Option::ignoreQV_);

    const bool isBaseMode = 
        options.is_set(Settings::Option::baseMode_) ? options.get(Settings::Option::baseMode_) : false;

    const bool isPulseMode = 
        options.is_set(Settings::Option::pulseMode_) ? options.get(Settings::Option::pulseMode_) : false;

    int modeCount = 0;
    if (isBaseMode)  modeCount++;
    if (isPulseMode) modeCount++;

    if (modeCount == 0)
        settings.mode = Settings::BaseMode;
    else if (modeCount == 1) 
        if (isBaseMode)  
            settings.mode = Settings::BaseMode;
        else 
            settings.mode = Settings::PulseMode;
    else
        settings.errors_.push_back("Unknown modes selected.");

    if (forcePulseMode) settings.mode = Settings::PulseMode;

    // BaseMap
    if (not options[Settings::Option::baseMap_].empty()) {
        settings.baseMap = options[Settings::Option::baseMap_];
        std::transform(settings.baseMap.begin(), settings.baseMap.end(), settings.baseMap.begin(), ::toupper);
        cout << settings.baseMap << endl;
        std::string _baseMap = settings.baseMap;
        std::sort(_baseMap.begin(), _baseMap.end());
        if (_baseMap != "ACGT") { settings.errors_.push_back("Bad basemap."); }
    }

    // input
    settings.inputBamFilenames = parser.args();
    if (settings.inputBamFilenames.size() == 1) {
        settings.polymeraseBamFilename = settings.inputBamFilenames[0];
        if (settings.polymeraseBamFilename.find("polymerase.bam") == std::string::npos)
            settings.errors_.push_back("missing input *.polymerase.bam.");
    } else if (settings.inputBamFilenames.size() == 2) {
        settings.subreadsBamFilename = settings.inputBamFilenames[0];
        settings.scrapsBamFilename   = settings.inputBamFilenames[1];
        if (settings.subreadsBamFilename.find("subreads.bam") == std::string::npos)
            settings.errors_.push_back("missing input *.subreads.bam.");
        if (settings.scrapsBamFilename.find("scraps.bam") == std::string::npos)
            settings.errors_.push_back("missing input *.scraps.bam.");
    } else {
        settings.errors_.push_back("missing input (polymerase.bam or subreads+scraps.bam.");
    }

    // output 
    settings.outputBaxPrefix = options[Settings::Option::output_];
    if (settings.outputBaxPrefix.empty()) { // if output prefix not set.
        if (not settings.subreadsBamFilename.empty()) {
            settings.outputBaxPrefix = internal::GetMovienameFromFilename(settings.subreadsBamFilename);
        } else if (not settings.polymeraseBamFilename.empty()) {
            settings.outputBaxPrefix = internal::GetMovienameFromFilename(settings.polymeraseBamFilename);
        }
    }

    if (settings.mode == Settings::BaseMode) 
        settings.outputBaxFilename = settings.outputBaxPrefix + ".bax.h5";
    else if (settings.mode == Settings::PulseMode)
        settings.outputBaxFilename = settings.outputBaxPrefix + ".plx.h5";
    settings.outputRgnFilename = settings.outputBaxPrefix + ".rgn.h5";

    // movie
    settings.movieName = internal::GetMovienameFromFilename(settings.outputBaxPrefix);

    if (options.is_set(Settings::Option::metadata_)) {
        // metadata.xml will be placed at upper directory of bax.h5
        settings.outputMetadataFilename = internal::DirectoryPath(settings.outputBaxPrefix) + "/../" + 
                                          settings.movieName + ".metadata.xml";
        settings.outputAnalysisDirname = internal::ParentDirectoryName(internal::DirectoryPath(settings.outputBaxPrefix));
    }

#ifdef DEBUG_SETTINGS
    string modeString = "Unknown";
    if (settings.mode == Settings::BaseMode)
        modeString = "base";
    else if (settings.mode == Settings::PulseMode)
        modeString = "pulse";

    cerr << "CommandLine: " << settings.program << " " << settings.args << endl
         << "Description: " << settings.description << endl
         << "Version    : " << settings.version << endl
         << "Mode       : " << modeString << endl
         << "BaseMap    : " << settings.baseMap << endl
         << "Movie name : " << settings.movieName << endl
         << "Input files: " << endl;
    if (not settings.subreadsBamFilename.empty())
         cerr << " subreads  : " << settings.subreadsBamFilename << endl;
    if (not settings.scrapsBamFilename.empty())
         cerr << " scraps    : " << settings.scrapsBamFilename << endl;
    if (not settings.polymeraseBamFilename.empty())
         cerr << " polymerase: " << settings.polymeraseBamFilename << endl;
    cerr << "Output h5  : " << settings.outputBaxFilename << endl
         << "Output xml : " << settings.outputMetadataFilename << endl;

#endif

    return settings;
}
