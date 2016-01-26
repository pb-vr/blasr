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

#include <unistd.h> // getcwd
#include <iostream>

#include "HDFBaxWriter.hpp"
#include "HDFPulseWriter.hpp"

#include "Bam2Bax.h"
#include "Bam2BaxConverter.h"
#include "IConverter.h"
#include <boost/scoped_ptr.hpp>

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
    boost::scoped_ptr<IConverter> converter;
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
