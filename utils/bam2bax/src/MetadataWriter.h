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

#ifndef _BAM2BAX_METADATA_WRITER_H_
#define _BAM2BAX_METADATA_WRITER_H_

#include <iostream>
#include <pbbam/ReadGroupInfo.h>

namespace internal{

const std::string DEFAULT_ANALYSIS_DIR = "Analysis_Results";

const std::string META_CONTENT = 
"<?xml version=\"1.0\" encoding=\"utf-8\"?><Metadata xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd=\"http://www.w3.org/2001/XMLSchema\" xmlns=\"http://pacificbiosciences.com/PAP/Metadata.xsd\"><InstCtrlVer>__BASECALLERVERSION__</InstCtrlVer><CellIndex>3</CellIndex><SetNumber>1</SetNumber><BindingKit><PartNumber>__BINDINGKIT__</PartNumber></BindingKit><SequencingKit><PartNumber>__SEQUENCINGKIT__</PartNumber></SequencingKit><Primary><Protocol>BasecallerV1</Protocol><ResultsFolder>__ANALYSISDIR__</ResultsFolder></Primary></Metadata>";

std::string Replace(const std::string & in_str,
                    const std::string & to_find,
                    const std::string & to_replace) {
    // Replace the first occurrence of to_find by to_replace.
    std::string ret = in_str;
    std::size_t pos = ret.find(to_find);
    if (pos != std::string::npos) {
        ret.replace(pos, to_find.size(), to_replace);
    }
    return ret;
}
} //namespace internal

class MetadataWriter {
public: 
    MetadataWriter(const std::string & filename, 
                   const PacBio::BAM::ReadGroupInfo & rg,
                   const std::string & analysisDir=internal::DEFAULT_ANALYSIS_DIR);

    MetadataWriter(const std::string & filename, 
                   const std::string & basecallerVersion,
                   const std::string & sequencingKit,
                   const std::string & bindingKit,
                   const std::string & analysisDir);

    ~MetadataWriter(void) {}
};

MetadataWriter::MetadataWriter(const std::string & filename, 
                               const PacBio::BAM::ReadGroupInfo & rg,
                               const std::string & analysisDir) {
    MetadataWriter(filename, 
                   rg.BasecallerVersion(),
                   rg.SequencingKit(),
                   rg.BindingKit(),
                   analysisDir);
}

MetadataWriter::MetadataWriter(const std::string & filename, 
                               const std::string & basecallerVersion,
                               const std::string & sequencingKit,
                               const std::string & bindingKit,
                               const std::string & analysisDir) {
    assert(analysisDir.find('/') == std::string::npos);
    ofstream ofile; 
    ofile.open(filename, std::ofstream::out);

    std::string to_print = internal::META_CONTENT;
    to_print = internal::Replace(to_print, "__BASECALLERVERSION__", basecallerVersion);
    to_print = internal::Replace(to_print, "__SEQUENCINGKIT__", sequencingKit);
    to_print = internal::Replace(to_print, "__BINDINGKIT__", bindingKit);
    to_print = internal::Replace(to_print, "__ANALYSISDIR__", analysisDir);

    ofile << to_print << endl;
    ofile.close();
}
#endif
