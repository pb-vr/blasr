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

// Author: Derek Barnett

#ifndef ICONVERTER_H
#define ICONVERTER_H

#include "Settings.h"
#include "SMRTSequence.hpp"
#include <pbbam/BamHeader.h>
#include <pbbam/BamWriter.h>
#include <map>
#include <string>
#include <vector>

namespace PacBio {
namespace BAM {

class BamRecordImpl;

} // namespace BAM
} // namespace PacBio

class IConverter
{
public:
    virtual ~IConverter(void);

public:
    virtual std::vector<std::string> Errors(void) const final;
    virtual bool Run(void) =0;

protected:
    IConverter(Settings& settings);

    virtual void AddErrorMessage(const std::string& e) final;

    virtual PacBio::BAM::BamHeader CreateHeader(const std::string& modeString) final;

    virtual std::string HeaderReadType(void) const =0;
    virtual std::string OutputFileSuffix(void) const =0;

protected:
    // common state
    Settings& settings_;
    std::vector<std::string> errors_;

    // run info for BamHeader creation
    std::string bindingKit_;
    std::string sequencingKit_;
    std::string basecallerVersion_;
    std::string frameRateHz_;
};

#endif // ICONVERTER_H
