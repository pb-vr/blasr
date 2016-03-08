// Author: Yuan Li

#include <unistd.h> // getcwd
#include <iostream>
#include <memory>

#include "HDFBaxWriter.hpp"
#include "HDFPulseWriter.hpp"

#include "Bam2Bax.h"
#include "Bam2BaxConverter.h"
#include "IConverter.h"

using namespace std;

namespace internal {
static inline
string CurrentWorkingDir(void)
{
    char result[FILENAME_MAX] = { };
    if (getcwd(result, FILENAME_MAX) == nullptr)
        return string();
    return string(result);
}

} // namespace internal

int Bam2Bax::Run(Settings& settings) {

    bool success = false;
    std::unique_ptr<IConverter> converter;
    if (settings.mode == Settings::BaseMode) {
        std::cout << "Converting BAM to bax.h5." << std::endl;
        converter.reset(new Bam2BaxConverter<HDFBaxWriter>(settings));
    }
    else if (settings.mode == Settings::PulseMode) {
        std::cout << "Converting BAM to plx.h5." << std::endl;
        converter.reset(new Bam2BaxConverter<HDFPulseWriter>(settings));
    }
    else { 
        cerr << "UNKNOWN mode." << settings.mode << endl;
        return EXIT_FAILURE;
    }

    if (converter->Run()) {
        success = true;
    }

    // return success/fail
    if (success)
        return EXIT_SUCCESS;
    else {
        for (const string& e : converter->Errors())
            cerr << "ERROR: " << e << endl;
        return EXIT_FAILURE;
    }
}
