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

#include "Bax2Bam.h"
#include "CcsConverter.h"
#include "HqRegionConverter.h"
#include "PolymeraseReadConverter.h"
#include "SubreadConverter.h"
#include <pbbam/DataSet.h>
#include <boost/algorithm/string.hpp>
#include <boost/scoped_ptr.hpp>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>

#include <unistd.h> // getcwd
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

static
bool WriteDatasetXmlOutput(const Settings& settings,
                           vector<string>* errors)
{
    using namespace PacBio::BAM;
    assert(errors);

    try {
        DataSet dataset(settings.datasetXmlFilename);
        assert(dataset.Type() == DataSet::HDF_SUBREAD);

        // change type
        dataset.Type(DataSet::SUBREAD);
        dataset.MetaType("PacBio.DataSet.SubreadSet");

        time_t currentTime = time(NULL);
        //const string& timestamp = CurrentTimestamp();
        dataset.CreatedAt(ToIso8601(currentTime));
        dataset.TimeStampedName(string{"pacbio_dataset_subreadset-"}+ToDataSetFormat(currentTime));

        // change files: remove BAX, add BAM
        std::vector<ExternalResource> toRemove;
        ExternalResources resources = dataset.ExternalResources();
        auto iter = resources.cbegin();
        auto end  = resources.cend();
        for (; iter != end; ++iter) {
            ExternalResource e = (*iter);
            boost::iterator_range<string::iterator> baxFound = boost::algorithm::ifind_first(e.MetaType(), "bax");
            if (!baxFound.empty()) 
                toRemove.push_back(e);
        }

        while(!toRemove.empty()) {
            auto e = toRemove.back();
            resources.Remove(e);
            toRemove.pop_back();
        }

        const string scheme = "file://";
        string fullOutputFilepath;

        // If the output filename starts with a slash, assume it's the path
        if (boost::starts_with(settings.outputBamFilename, "/"))
        {
            fullOutputFilepath = settings.outputBamFilename;
        }
        else // otherwise build the path from the CWD
        { 
            fullOutputFilepath = CurrentWorkingDir();
            if (!fullOutputFilepath.empty())
                fullOutputFilepath.append(1, '/');
            fullOutputFilepath.append(settings.outputBamFilename);
        }

        // Combine the scheme and filepath and store in the dataset
        fullOutputFilepath = scheme + fullOutputFilepath;
        resources.Add(ExternalResource("SubreadFile.SubreadBamFile",fullOutputFilepath));
        dataset.ExternalResources(resources);

        // save to file 
        string xmlFn = settings.outputXmlFilename; // try user-provided explicit filename first
        if (xmlFn.empty())
            xmlFn = settings.outputBamPrefix + ".dataset.xml"; // prefix set w/ moviename elsewhere if not user-provided
        dataset.Save(xmlFn);
        return true;

    } catch (std::exception&) {
        errors->push_back("could not create output XML");
        return false;
    }
}

} // namespace internal

int Bax2Bam::Run(Settings& settings) {

    // init conversion mode
    boost::scoped_ptr<IConverter> converter;
    switch (settings.mode) {
        case Settings::HQRegionMode   : converter.reset(new HqRegionConverter(settings)); break;
        case Settings::PolymeraseMode : converter.reset(new PolymeraseReadConverter(settings)); break;
        case Settings::SubreadMode    : converter.reset(new SubreadConverter(settings)); break;
        case Settings::CCSMode        : converter.reset(new CcsConverter(settings)); break;
        default :
            cerr << "ERROR: unknown mode selected" << endl;
            return EXIT_FAILURE;
    }

    // run conversion
    bool success = false;
    vector<string> xmlErrors;
    if (converter->Run()) {
        success = true;

        // if given dataset XML as input, attempt write dataset XML output
        if (!settings.datasetXmlFilename.empty()) {
            if (!internal::WriteDatasetXmlOutput(settings, &xmlErrors))
                success = false;
        }
    }

    // return success/fail
    if (success)
        return EXIT_SUCCESS;
    else {
        for (const string& e : converter->Errors())
            cerr << "ERROR: " << e << endl;
        for (const string& e : xmlErrors)
            cerr << "ERROR: " << e << endl;
        return EXIT_FAILURE;
    }
}
