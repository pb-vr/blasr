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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <string>
#include <vector>

namespace optparse { class OptionParser; }

class Settings
{
public:
    enum Mode { BaseMode  // BAM to BAX.H5
              , PulseMode // BAM to PLS.H5
              }; 

public:
    Settings(void);
    static Settings FromCommandLine(optparse::OptionParser& parser,
                                    int argc,
                                    char* argv[],
                                    bool forcePulseMode=false);
    struct Option {
        static const char* input_;
        static const char* output_;
        static const char* metadata_;
        static const char* baseMode_;
        static const char* pulseMode_;
        static const char* ignoreQV_;
        static const char* baseMap_;
    };

    // default option value
    struct OptionValue {
        static const char* baseMap_;
    };

public:
    // input
    std::vector<std::string> inputBamFilenames;

    std::string subreadsBamFilename;
    std::string scrapsBamFilename;
    std::string polymeraseBamFilename;

    //output
    std::string outputBaxPrefix;
    std::string outputBaxFilename;
    std::string outputRgnFilename;

    std::string outputMetadataFilename;
    std::string outputAnalysisDirname;

    // program info
    std::string program;
    std::string args;
    std::string version;
    std::string description;

    // generated
    std::string movieName;

    // mode
    Mode mode;

    bool ignoreQV;

    // base map
    std::string baseMap;

    // command line parsing
    std::vector<std::string> errors_;
};

#endif // SETTINGS_H
