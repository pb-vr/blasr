// Author: Derek Barnett

#include "TestData.h"
#include "TestUtils.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <cstdlib>

using namespace std;

void RemoveFiles(const vector<string>& filenames)
{
    for (auto fn : filenames)
        remove(fn.c_str());
}

void RemoveFile(const string& filename)
{
    vector<string> filenames;
    filenames.push_back(filename);
    RemoveFiles(filenames);
}

int RunBam2Bax(const vector<string>& bamFilenames,
               const string& outputType,
               const string& additionalArgs)
{
    string convertArgs;
    convertArgs += outputType;
    if (!additionalArgs.empty()) {
        convertArgs += string(" ");
        convertArgs += additionalArgs;
    }
    for (auto fn : bamFilenames) {
        convertArgs += string(" ");
        convertArgs += fn;
    }

    const string& convertCommandLine = tests::Bam2Bax_Exe + string(" ") + convertArgs;
    return system(convertCommandLine.c_str());
}
