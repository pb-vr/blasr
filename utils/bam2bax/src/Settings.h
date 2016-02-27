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
