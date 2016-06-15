// Author: Derek Barnett

#include "Bax2Bam.h"
#include "CcsConverter.h"
#include "HqRegionConverter.h"
#include "PolymeraseReadConverter.h"
#include "SubreadConverter.h"
#include <pbbam/DataSet.h>
#include <pbbam/PbiRawData.h>
#include <boost/algorithm/string.hpp>
#include <memory>
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

        // determine output details based on mode
        // initialize with SUBREAD data (most common)
        DataSet::TypeEnum outputDataSetType;
        string outputDataSetMetaType;
        string outputTimestampPrefix;
        string outputBamFileType;
        string outputScrapsFileType;
        string outputXmlSuffix;

        switch(settings.mode)
        {
            case Settings::SubreadMode :
            {
                outputDataSetType = DataSet::SUBREAD;
                outputDataSetMetaType = "PacBio.DataSet.SubreadSet";
                outputTimestampPrefix = "pacbio_dataset_subreadset-";
                outputBamFileType = "PacBio.SubreadFile.SubreadBamFile";
                outputScrapsFileType = "PacBio.SubreadFile.ScrapsBamFile";
                outputXmlSuffix = ".subreadset.xml";
                break;
            }

            case Settings::CCSMode :
            {
                outputDataSetType = DataSet::CONSENSUS_READ;
                outputDataSetMetaType = "PacBio.DataSet.ConsensusReadSet";
                outputTimestampPrefix = "pacbio_dataset_consensusreadset-";
                outputBamFileType = "PacBio.ConsensusReadFile.ConsensusReadBamFile";
                outputScrapsFileType = "";
                outputXmlSuffix = ".consensusreadset.xml";
                break;

            }
            case Settings::HQRegionMode :
            {
                outputDataSetType = DataSet::SUBREAD;
                outputDataSetMetaType = "PacBio.DataSet.SubreadSet";
                outputTimestampPrefix = "pacbio_dataset_subreadset-";
                outputBamFileType = "PacBio.SubreadFile.HqRegionBamFile";
                outputScrapsFileType = "PacBio.SubreadFile.HqScrapsBamFile";;
                outputXmlSuffix = ".subreadset.xml";
                break;
            }
            case Settings::PolymeraseMode :
            {
                outputDataSetType = DataSet::SUBREAD;
                outputDataSetMetaType = "PacBio.DataSet.SubreadSet";
                outputTimestampPrefix = "pacbio_dataset_subreadset-";
                outputBamFileType = "PacBio.SubreadFile.PolymeraseBamFile";
                outputScrapsFileType = "PacBio.SubreadFile.PolymeraseScrapsBamFile";
                outputXmlSuffix = ".subreadset.xml";
                break;
            }

            default:
                assert(false); // should already be checked upstream
                throw std::runtime_error("unknown mode selected");
        }

        DataSet dataset(settings.datasetXmlFilename);
        assert(dataset.Type() == DataSet::HDF_SUBREAD);

        // change type
        dataset.Type(outputDataSetType);
        dataset.MetaType(outputDataSetMetaType);

        time_t currentTime = time(NULL);
        //const string& timestamp = CurrentTimestamp();
        dataset.CreatedAt(ToIso8601(currentTime));
        dataset.TimeStampedName(outputTimestampPrefix+ToDataSetFormat(currentTime));

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
        string mainBamFilepath;

        // If the output filename starts with a slash, assume it's the path
        if (boost::starts_with(settings.outputBamFilename, "/"))
        {
            mainBamFilepath = settings.outputBamFilename;
        }
        else // otherwise build the path from the CWD
        { 
            mainBamFilepath = CurrentWorkingDir();
            if (!mainBamFilepath.empty())
                mainBamFilepath.append(1, '/');
            mainBamFilepath.append(settings.outputBamFilename);
        }

        // Combine the scheme and filepath and store in the dataset
        mainBamFilepath = scheme + mainBamFilepath;
        ExternalResource mainBam{ outputBamFileType, mainBamFilepath };
        FileIndex mainPbi{ "PacBio.Index.PacBioIndex", mainBamFilepath + ".pbi" };
        mainBam.FileIndices().Add(mainPbi);

        // maybe add scraps BAM (& PBI)
        if (!settings.scrapsBamFilename.empty()) {

            string scrapsBamFilepath;

            // If the output filename starts with a slash, assume it's the path
            if (boost::starts_with(settings.scrapsBamFilename, "/"))
            {
                scrapsBamFilepath = settings.scrapsBamFilename;
            }
            else // otherwise build the path from the CWD
            {
                scrapsBamFilepath = CurrentWorkingDir();
                if (!scrapsBamFilepath.empty())
                    scrapsBamFilepath.append(1, '/');
                scrapsBamFilepath.append(settings.scrapsBamFilename);
            }

            ExternalResource scrapsBam{ outputScrapsFileType, scrapsBamFilepath };
            FileIndex scrapsPbi{ "PacBio.Index.PacBioIndex", scrapsBamFilepath + ".pbi" };
            scrapsBam.FileIndices().Add(scrapsPbi);
            mainBam.ExternalResources().Add(scrapsBam);
        }

        // add resources to output dataset
        resources.Add(mainBam);
        dataset.ExternalResources(resources);

        // update TotalLength & NumRecords
        const BamFile subreadFile{ settings.outputBamFilename };
        const string subreadPbiFn = subreadFile.PacBioIndexFilename();
        const PbiRawData subreadsIndex{ subreadPbiFn };
        const PbiRawBasicData& subreadData = subreadsIndex.BasicData();

        uint64_t totalLength = 0;
        uint32_t numRecords = subreadsIndex.NumReads();
        for (uint32_t i = 0; i < numRecords; ++i) {
            const auto subreadLength = subreadData.qEnd_.at(i) - subreadData.qStart_.at(i);
            totalLength += subreadLength;
        }

        DataSetMetadata metadata = dataset.Metadata();
        metadata.TotalLength(std::to_string(totalLength));
        metadata.NumRecords(std::to_string(numRecords));
        dataset.Metadata(metadata);

        // save to file 
        string xmlFn = settings.outputXmlFilename; // try user-provided explicit filename first
        if (xmlFn.empty())
            xmlFn = settings.outputBamPrefix + outputXmlSuffix; // prefix set w/ moviename elsewhere if not user-provided
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
    std::unique_ptr<IConverter> converter;
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
