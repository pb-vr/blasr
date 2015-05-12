#include "utils/FileOfFileNames.hpp"
#include "datastructures/alignmentset/SAMSupplementalQVList.hpp"
#include "format/SAMHeaderPrinter.hpp"
#include "format/BAMPrinter.hpp"
#include "pbbam/BamWriter.h"
#include "CommandLineParser.hpp"

using namespace PacBio::BAM;
using namespace std;

string DISCLAIM = 
"THIS TOOL IS CREATED FOR DEVELOPERS USE ONLY AND IT MAY OR MAY NOT "
"BREAK AT ANY TIME. USE AT YOUR OWN RISK.";

string GetVersion(void) {
    return "1.0";
}

void CCSReadToBamRecord(CCSSequence & ccsRead, BamRecord & bamRecord, SupplementalQVList & samQVList) {
//m140905_042212_sidney_c100564852550000001823085912221377_s1_X0/1650/1920_2155   4   *   0   255 *   *   0   0
    bamRecord.Impl().Name(ccsRead.GetTitle());
    bamRecord.Impl().Flag(static_cast<uint32_t>(4));
    string seqString;
    seqString.assign((char*)ccsRead.seq, ccsRead.length);

    bamRecord.Impl().SetSequenceAndQualities(seqString, ccsRead.qual.ToString());
    // bamRecord.Impl().CigarData(Cigar::FromStdString("*"));
    bamRecord.Impl().Bin(0);
    bamRecord.Impl().InsertSize(0);
    bamRecord.Impl().MatePosition(static_cast<PacBio::BAM::Position>(-1));
    bamRecord.Impl().MateReferenceId(static_cast<int32_t>(-1));
    bamRecord.Impl().Position(static_cast<PacBio::BAM::Position>(-1));
    bamRecord.Impl().ReferenceId(static_cast<int32_t>(-1));
    TagCollection tags;
    tags["RG"] = ccsRead.GetReadGroupId();
    tags["np"] = ccsRead.numPasses;
    tags["zm"] = ccsRead.zmwData.holeNumber;
    tags["qs"] = 0;
    tags["qe"] = ccsRead.length;

    samQVList.FormatQVOptionalFields(ccsRead);
    // Add QVs to BamRecordImpl.
    string insertionQVs, deletionQVs, substitutionQVs, mergeQVs, substitutionTags, deletionTags;
    if (ccsRead.GetQVs("InsertionQV", insertionQVs)) {
        tags["iq"] = insertionQVs;
    }
    if (ccsRead.GetQVs("DeletionQV", deletionQVs)) {
        tags["dq"] = deletionQVs;
    }
    if (ccsRead.GetQVs("SubstitutionQV", substitutionQVs)) {
        tags["sq"] = substitutionQVs;
    }
    if (ccsRead.GetQVs("MergeQV", mergeQVs)) {
        tags["mq"] = mergeQVs;
    }
    // substitutionTag is not included by default
    if (ccsRead.GetQVs("DeletionTag", deletionTags)) {
        tags["dt"] = deletionTags;
    }
    bamRecord.Impl().Tags(tags);
}

int main(int argc, char* argv[]) {
    string progName = "ccsh5tobam";
    CommandLineParser clp;
    clp.SetHelp("Convert ccs.h5 to bam.\n" + DISCLAIM);
    clp.SetConciseHelp("ccsh5tobam ccs.h5|fofn out.bam\n" + DISCLAIM);
    clp.SetProgramName(progName);
    clp.SetVersion(GetVersion());
    string fofn, bamOutName;
    clp.RegisterStringOption("in.ccs.h5", &fofn, "Input ccs.h5|fofn file.", true);
    clp.RegisterStringOption("out.bam", &bamOutName, "Output bam file.", true);
    clp.RegisterPreviousFlagsAsHidden();
    clp.ParseCommandLine(argc, argv);

    //cerr << "[INFO] " << GetTimestamp() << " [" << progName << "] started."  << endl;

    vector<string> ccsFileNames;
    FileOfFileNames::StoreFileOrFileList(fofn, ccsFileNames);

    string so = "UNKNOWN"; // sorting order;
    string version = GetVersion(); 
    string commandLineString;
    clp.CommandLineToString(argc, argv, commandLineString);

    SupplementalQVList samQVList;
    samQVList.SetDefaultQV();
    SequenceIndexDatabase<FASTASequence> seqdb;
    SAMHeaderPrinter shp(so, seqdb,
            ccsFileNames, ReadType::ReadTypeEnum::CCS,
            samQVList, "ccsh52bam", version,
            commandLineString);
    string headerString = shp.ToString();// SAM/BAM header

    BamHeader header = BamHeader(headerString);
    // Both file name and SAMHeader are required in order to create a BamWriter.
    BamWriter * bamWriterPtr = new BamWriter(bamOutName, header);

    for (string ccsFileName: ccsFileNames) {
        ReaderAgglomerate reader;
        reader.SetReadFileName(ccsFileName);
        reader.SetReadType(ReadType::ReadTypeEnum::CCS);

        // Initialize using already set file names.
        int initReturnValue = reader.Initialize();
        if (initReturnValue <= 0) {
            cerr << "WARNING! Could not open file " << ccsFileName << endl;
            continue;
        }

        // Check whether use ccs only.
        assert (reader.GetFileType() == HDFCCSONLY);
        int randint = 0;
        CCSSequence ccsRead;
        while(reader.GetNext(ccsRead, randint) != 0) {
            if (ccsRead.length > 0) {
                BamRecord bamRecord;
                CCSReadToBamRecord(ccsRead, bamRecord, samQVList);
                bamWriterPtr->Write(bamRecord);
            }
        }
    }

    try {
        bamWriterPtr->TryFlush();
        delete bamWriterPtr;
        bamWriterPtr = NULL;
    } catch (std::exception e) {
        cout << "Error, could not flush bam records to bam file." << endl;
        exit(1);
    }

    //cerr << "[INFO] " << GetTimestamp() << " [" << progName << "] ended."  << endl;
    return 0;
}
