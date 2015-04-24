#include <vector>
#include <string>
#include "NucConversion.hpp"
#include "FASTASequence.hpp"
#include "FASTAReader.hpp"
#include "suffixarray/SuffixArray.hpp"
#include "suffixarray/SuffixArrayTypes.hpp"
#include "suffixarray/ssort.hpp"
#include "algorithms/sorting/qsufsort.hpp"


void PrintUsage() {
    cout << "samodify changes word size of input suffix array." << endl;
	cout << "Usage: samodify in.sa genome.fasta out.sa [-blt p]" << endl;
	cout << "       -blt p  Build a lookup table on prefixes of length 'p' " << endl;
}

int main(int argc, char* argv[]) {

    if (argc < 4) {
        PrintUsage();
        exit(1);
    }
    int argi = 1;
    string saInFile = argv[argi++];
    string genomeFileName = argv[argi++];
    string saOutFile = argv[argi++];
    vector<string> inFiles;

    int doBLT = 0;
    int doBLCP = 0;
    int bltPrefixLength = 0;
    int lcpLength = 0;
    int parsingOptions = 0;

    while (argi < argc) {
        if (strcmp(argv[argi], "-blt") == 0) {
            doBLT = 1;
            bltPrefixLength = atoi(argv[++argi]);
        }
        else if (strcmp(argv[argi], "-blcp") == 0) {
            doBLCP = 1;
            lcpLength = atoi(argv[++argi]);
        }
        else {
            PrintUsage();
            cout << "Bad option: " << argv[argi] << endl;
            exit(1);
        }
        ++argi;
    }

    //
    // Read the suffix array to modify.
    //

    DNASuffixArray  sa;
    sa.Read(saInFile);

    FASTAReader reader;
    reader.Initialize(genomeFileName);
    FASTASequence seq;
    reader.ReadAllSequencesIntoOne(seq);


    if (doBLT) {
        sa.BuildLookupTable(seq.seq, seq.length, bltPrefixLength);
    }

    if (doBLCP) {
        cout << "LCP Table not yet implemented." << endl;
    }

    sa.Write(saOutFile);

}

