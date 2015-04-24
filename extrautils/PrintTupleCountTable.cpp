#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "utils.hpp"
#include "FASTASequence.hpp"
#include "FASTAReader.hpp"
#include "CommandLineParser.hpp"
#include "tuples/DNATuple.hpp"
#include "tuples/CompressedDNATuple.hpp"
#include "tuples/TupleMetrics.hpp"
#include "tuples/TupleCountTable.hpp"


#ifdef COMPRESSED
typedef TupleCountTable<FASTASequence, CompressedDNATuple<FASTASequence> > CountTable;
#else
typedef TupleCountTable<FASTASequence, DNATuple> CountTable;
#endif

int main(int argc, char* argv[]) {
    CommandLineParser clp;
    string tableFileName;
    vector<string> sequenceFiles;
    TupleMetrics tm;
    int tupleSize = 8;
    clp.SetProgramName("printTupleCountTable");
    clp.SetProgramSummary("Count the number of occurrences of every k-mer in a file.");
    clp.RegisterStringOption("table", &tableFileName, "Output table name.", true);
    clp.RegisterIntOption("wordsize", &tupleSize, "Size of words to count", 
                          CommandLineParser::NonNegativeInteger, false);
    clp.RegisterStringListOption("reads", &sequenceFiles, "All sequences.", false);
    clp.RegisterPreviousFlagsAsHidden();
    vector<string> opts;

    if (argc == 2) {
        string fastaFileName = argv[1];
        sequenceFiles.push_back(fastaFileName);
        tableFileName = fastaFileName + ".ctab";
    }
    else {
        clp.ParseCommandLine(argc, argv, opts);
    }

    tm.tupleSize = tupleSize;
    tm.InitializeMask();
    ofstream tableOut;
    CrucialOpen(tableFileName, tableOut, std::ios::out| std::ios::binary);
    CountTable table;
    table.InitCountTable(tm);
    int i;
    FASTASequence seq;
    for (i = 0; i < sequenceFiles.size(); i++ ){ 
        FASTAReader reader;
        reader.Init(sequenceFiles[i]);
        while (reader.GetNext(seq)) {
            seq.ToUpper();
            table.AddSequenceTupleCountsLR(seq);
        }
    }	
    table.Write(tableOut);

    return 0;
}
