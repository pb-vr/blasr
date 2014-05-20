/*
 * =====================================================================================
 *
 *       Filename:  SAMFilter.cpp
 *
 *    Description:  Filter SAM Hits according to 
 *                  filteration criteria
 *                     minPctSimilarity, minAccuracy, 
 *                     minLength, holeNumbers
 *                  and multiple-hit policy
 *                     random    : a random hit
 *                     all       : all hits
 *                     allbest   : all hits with the best score
 *                     randombest: a random hit selected from all the hits 
 *                                 that have the best score
 *
 *        Version:  1.0
 *        Created:  03/19/2013 01:19:43 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * =====================================================================================
 */
#include <iostream>
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "CommandLineParser.h"
#include "utils/ChangeListID.h"
#include "utils/TimeUtils.h"
#include "utils/RangeUtils.h"
#include "utils/SMRTReadUtils.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "algorithms/alignment/AlignmentUtils.h"
#include "algorithms/alignment/StringToScoreMatrix.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "algorithms/alignment/printers/SAMPrinter.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/FilterCriteria.h"
#include "datastructures/metagenome/TitleTable.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
#include "data/gff/GFFFile.h"

//#define USE_GOOGLE_PROFILER
#ifdef USE_GOOGLE_PROFILER
#include "gperftools/profiler.h"
#endif

char VERSION[] = "v0.1.0";
char PERFORCE_VERSION_STRING[] = "$Change: 134995 $";
// By default negative score is better.
SCORESIGN scoreSign = NEG;

// Score functions can be: 
//   ALIGNERSCORE: aligner's score as indicated by 'as' flag
//   EDITDIST    : edit distance between query and target
//   BLASESCORE  : blasr's default score
//   USERSCORE   : score computed from user specifed score matrix, and
//                 insertion & deletion scores
enum SCOREFUNC {ALIGNERSCORE, EDITDIST, BLASRSCORE, USERSCORE};

enum HITPOLICY {RANDOM, ALL, ALLBEST, RANDOMBEST, LEFTMOST};

// Return a multiple hit policy.
HITPOLICY setHitPolicy(const string & hitPolicyStr) {
    if (hitPolicyStr == "random")          return RANDOM;
    else if (hitPolicyStr == "all")        return ALL;
    else if (hitPolicyStr == "allbest")    return ALLBEST;
    else if (hitPolicyStr == "randombest") return RANDOMBEST;
    else if (hitPolicyStr == "leftmost")   return LEFTMOST;
    else {
        cout <<"ERROR, the specified multiple hit policy " 
             << hitPolicyStr <<" is not supported." << endl;
        exit(1);
    }
}

// Return a score function for computing alignment scores.
SCOREFUNC setScoreFunction(const string & scoreFuncStr) {
    if (scoreFuncStr == "alignerscore")    return ALIGNERSCORE;
    else if (scoreFuncStr == "editdist")   return EDITDIST;
    else if (scoreFuncStr == "blasrscore") return BLASRSCORE;
    else if (scoreFuncStr == "userscore")  return USERSCORE;
    else {
        cout <<"ERROR, the specified score function " 
        << scoreFuncStr <<" is not supported." << endl;
        exit(1);
    }
}

// Compare SAMAlignment objects by qName, score and 
// target positions.
bool byQNameScoreTStart(const SAMAlignment & a, 
                        const SAMAlignment & b) {
    if (a.qName == b.qName) {
        if (a.score == b.score) 
            return a.pos < b.pos;
        return Score(a.score, scoreSign).WorseThan(b.score);
    }
    return (a.qName < b.qName);
}

// Compare SAMAlignment objects by rName and qName 
bool byRNameQName(const SAMAlignment & a, 
                  const SAMAlignment & b) {
    if (a.rName == b.rName) {
        return a.qName < b.qName;
    }
    return (a.rName < b.rName);
}

// Get the next group of SAM alignments that have the same qName from
// allSAMAlignments[groupBegin ... groupEnd)
// Note that allSAMAlignments is already sorted by qName, score and tPos.
void GetNextSAMAlignmentGroup(vector<SAMAlignment> & allSAMAlignments, 
                              unsigned int groupBegin, 
                              unsigned int & groupEnd) {
    assert(groupBegin < allSAMAlignments.size());
    groupEnd = groupBegin + 1;
    string queryName = allSAMAlignments[groupBegin].qName;
    while(groupEnd < allSAMAlignments.size()) {
        if (allSAMAlignments[groupEnd].qName == queryName) 
            groupEnd ++;
        else break;
    }
}

// Get the best SAM alignments whose alignment score are the best. 
// Assume that alignments in allSAMAlignments[groupBegin, groupEnd)
// all have the same queryName and are sorted by score and tPos 
// asscendingly: worst, ...., best
void GetBestSAMAlignmentsInGroup(vector<SAMAlignment> & allSAMAlignments, 
                                 const unsigned int & groupBegin, 
                                 const unsigned int & groupEnd, 
                                 unsigned int & bestBegin,
                                 unsigned int & bestEnd) {
    assert(groupEnd  <= allSAMAlignments.size() and
           groupBegin < groupEnd);

    bestEnd = groupEnd;
    bestBegin = groupEnd - 1;
    int groupBestScore = allSAMAlignments[bestBegin].score;
    string queryName = allSAMAlignments[bestBegin].qName;
    while (bestBegin >= groupBegin and bestBegin < groupEnd) {
        assert(allSAMAlignments[bestBegin].qName == queryName);
        if (allSAMAlignments[bestBegin].score == groupBestScore)
            bestBegin -= 1;
        else break;
    }
    bestBegin += 1;
}

// Apply hit policy to a group of SAM alignments and return indices
// of the selected alignments.
vector<unsigned int> ApplyHitPolicy(HITPOLICY hitPolicy, 
                                    vector<SAMAlignment> & allSAMAlignments, 
                                    const unsigned int & groupBegin, 
                                    const unsigned int & groupEnd) {
    vector<unsigned int> hitIndices;
    if (hitPolicy == ALL) {
        for(unsigned int i = groupBegin; i < groupEnd; i++){
            hitIndices.push_back(i);
        }
    } else if (hitPolicy == RANDOM) {
        hitIndices.push_back(rand()%(groupEnd - groupBegin) + groupBegin);
    } else {
        unsigned int bestBegin, bestEnd;
        GetBestSAMAlignmentsInGroup(allSAMAlignments, groupBegin, groupEnd,
                                    bestBegin, bestEnd);
        if (hitPolicy == ALLBEST) {
            for(unsigned int i = bestBegin; i < bestEnd; i++){
                hitIndices.push_back(i);
            }
        } else if (hitPolicy == RANDOMBEST) {
            hitIndices.push_back(rand()%(bestEnd-bestBegin) + bestBegin);
        } else if (hitPolicy == LEFTMOST) {
            hitIndices.push_back(bestBegin);
        } else {
            assert(false);
        }
    }
    return hitIndices;
}

// Convert references[...].title in reference.fasta to their corresponding
// indices in the title table.
void ConvertTitlesToTitleTableIndices(vector<FASTASequence> & references,
        string & titleTableName) {
    TitleTable tt;
    tt.Read(titleTableName);
    for(int i = 0; i < references.size(); i++) {
        string title = references[i].GetTitle();
        int idx = -1;
        if (tt.Lookup(title, idx)) {
            stringstream ss;
            ss << idx; 
            references[i].CopyTitle(ss.str());
        } else {
            cout << "ERROR, reference " << title << " does not exist "
                 << " in the title table " << titleTableName << ". The "
                 << "reference fasta and the title table do not match."
                 << endl;
            exit(1);
        }
    }
    tt.Free();
}

// Return true if the alignment can only map to an adapter specified 
// in the adapter GFF file. 
// A sample record in adapter GFF file:
// ref000001   .   adapter 10955   10999   0.00    +   .   xxxx
// ref000001   .   adapter 32886   32930   0.00    +   .   xxxx 
// Note that the first field (e.g., 'ref000001') is id of sequence  
// in a reference repository, not sequence name, so we need to 
// reconstruct the mapping between sequence id and sequence name.
bool CheckAdapterOnly(GFFFile & adapterGffFile, //Adapter gff file
    AlignmentCandidate<> & alignment, // An alignment
    map<string, int> & refNameToIndex) { 
    // Map target sequence name to its index in reference repository.
    if (refNameToIndex.find(alignment.tName) == refNameToIndex.end()) {
        // This should not happen ...
        cout << "ERROR, could not find alignment target name "
             << alignment.tName << " in the reference file." << endl;
        exit(1);
    }
    int refNameIndex = refNameToIndex[alignment.tName];
    char buf [16];
    sprintf(buf, "ref%06d", refNameIndex + 1);
    // Reconstruct ref id in the format "ref00000?".
    string refNameId(buf);
    int FUZZY_OVERLAP = 20;
    for(int eindex = 0; eindex < adapterGffFile.entries.size();
            eindex++) { 
        GFFEntry & entry = adapterGffFile.entries[eindex];
        // Convert each GFF record from 1-based inclusive to 
        // 0-based exclusive.
        if (entry.type == "adapter" and 
            (entry.name == alignment.tName or
             entry.name == refNameId)) {
            UInt estart = entry.start - 1;
            UInt eend = entry.end;
            if (entry.strand == '-') {
                UInt tmp = estart;
                estart = alignment.tLength - 1 - eend;
                eend = alignment.tLength - 1 - tmp;
            }
            if (not (eend < alignment.GenomicTBegin() or
                 estart > alignment.GenomicTEnd())) {
                int lengthUnion = max(eend, alignment.GenomicTEnd()) -
                                  min(estart, alignment.GenomicTBegin());
                if (lengthUnion < eend - estart + FUZZY_OVERLAP) {
                    return true;
                }
            }
        }
    }
    return false;
}

int main(int argc, char* argv[]) {
#ifdef USE_GOOGLE_PROFILER
    char *profileFileName = getenv("CPUPROFILE");
    if (profileFileName != NULL) {
      ProfilerStart(profileFileName);
    }
    else {
      ProfilerStart("google_profile.txt");
    }
#endif
    string samFileName, refFileName, outFileName;
    int  scoreCutoff    = INF_INT;
    int  scoreSignInt   = -1;
    string scoreFuncStr = "alignerscore";
    string scoreMatrixStr= "";
    SCOREFUNC scoreFunc;

    int insScore = 5; // Insertion penalty
    int delScore = 5; // Deletion penalty

    int  seed           = 1; 
    bool isSorted       = false; // Whether input sam file is sorted.
    int  verbosity      = 0;

    string holeNumberStr;
    Ranges holeNumberRanges;

    FilterCriteria filterCriteria;
    string hitPolicyStr = "randombest";
    HITPOLICY hitPolicy = RANDOMBEST; 

    bool parseSmrtTitle = false;
    string titleTableName = "";
    string adapterGffFileName = "";

    CommandLineParser clp;
    clp.RegisterStringOption("file.sam", &samFileName,
                             "Input SAM file.");
    clp.RegisterStringOption("reference.fasta", &refFileName,
                             "Reference used to generate reads.");
    clp.RegisterStringOption("out.sam", &outFileName,
                             "Output SAM file.");

    clp.RegisterPreviousFlagsAsHidden();

    stringstream help;
    help << "(" << filterCriteria.minAccuracy 
         << ") Minimum percentage accuracy to reference.";
    clp.RegisterFloatOption("minAccuracy", &filterCriteria.minAccuracy,
            help.str(), CommandLineParser::PositiveFloat);

    help.str(string());
    help << "(" << filterCriteria.minPctSimilarity
         << ") Minimum percentage similarity to reference.";
    clp.RegisterFloatOption("minPctSimilarity", &filterCriteria.minPctSimilarity,
            help.str(), CommandLineParser::PositiveFloat);

    help.str(string());
    help << "(" << filterCriteria.minLength
         << ") Minimum aligned read length to output a hit."; 
    clp.RegisterIntOption("minLength", &filterCriteria.minLength,
            help.str(), CommandLineParser::PositiveInteger);
    clp.RegisterStringOption("hitPolicy", &hitPolicyStr,
            "(randombest) Specify a policy to treat multiple hits from "
            "[random, all, allbest, randombest]\n"
            "  random  : selects a random hit.\n"
            "  all     : selects all hits.\n"
            "  allbest : selects all the best alignment score hits.\n"
            "  randombest: selects a random hit from all best alignment score hits.\n"
            "  leftmost: selects a hit which has the best alignment score and \n" 
            "            has the smallest mapping coordinate in any reference.");
    clp.RegisterStringOption("scoreFunction", &scoreFuncStr,
            "(alignerscore) Specify an alignment score function from "
            "[alignerscore, editdist, blasrscore, userscore]\n" // affine
            "  alignerscore : aligner's score in SAM tag 'as'.\n"
            "  editdist     : edit distance between read and reference.\n"
            "  blasrscore   : blasr's default score function.\n"
            "  userscore    : score computed using a user-defined scoring\n"
            "                 matrix (specified by -scoreMatrix) and \n"
            "                 insertion & deletion scores (specified by \n"
            "                 -insertion and -deletion respectively).");
    clp.RegisterStringOption("scoreMatrix", &scoreMatrixStr,
            "Specify a user-defined score matrix for scoring reads."
            "The matrix is in the format\n"
            "    ACGTN\n"
            "  A abcde\n"
            "  C fghij\n"
            "  G klmno\n"
            "  T pqrst\n"
            "  N uvwxy\n" 
            ". The values a...y should be input as a quoted space "
            "seperated string: \"a b c ... y\". Lower scores are better, "
            "so matches should be less than mismatches e.g. a,g,m,s = -5 "
            "(match), mismatch = 6. "); 
    clp.RegisterIntOption("deletion", &delScore, 
            "Specify a user-defined deletion score.",
            CommandLineParser::Integer);
    clp.RegisterIntOption("insertion", &insScore,
            "Specify a user-defined insertion score.",
            CommandLineParser::Integer); 
    clp.RegisterIntOption("scoreCutoff", &scoreCutoff,
            "Score cut off defines the worst score to output a hit.", 
            CommandLineParser::Integer);
    clp.RegisterIntOption("scoreSign", &scoreSignInt,
            "(-1) Sepcifiy the alignment score sign."
            "  -1: negative scores are better than positive ones."
            "   1: positive scores are better than negative ones.",
            CommandLineParser::Integer);
    clp.RegisterIntOption ("seed", &seed,
            "(1)  Seed for random number generator."
            "   If seed is 0, then the current time will be used as seed.",
            CommandLineParser::Integer);
    clp.RegisterStringOption("holeNumbers", &holeNumberStr,
            "A string of comma-delimited hole number ranges to output hits, "
            "such as '1,2,10-12'. "
            "This requires hit titles to be in SMRT read title format.");
    clp.RegisterFlagOption("smrtTitle", &parseSmrtTitle,
            "Use this option when filtering alignments generated by "
            "programs other than blasr, e.g. bwa-sw or gmap. "
            "  Parse read coordinates from the SMRT read title. " 
            "The title is in the format /name/hole/coordinates, where"
            " coordinates are in the format \\d+_\\d+, and represent "
            "the interval of the read that was aligned.");
    /* This experimental option can be useful for metagenomics, in which case
     * there are hundreds of sequences in the target, of which many titles are
     * long and may contain white spaces (e.g., ' ', '\t'). 
     * In order to save disc space and avoid the (possibly) none unique mapping
     * between full and short reference names, one may call blasr with 
     * -titleTable option to represent all target sequences in the output
     * by their indices in the title table.*/
    clp.RegisterStringOption("titleTable", &titleTableName,
            "Use this experimental option when filtering alignments generated by "
            "blasr with -titleTable titleTableName, in which case "
            "reference titles in SAM are represented by their "
            "indices (e.g., 0, 1, 2, ...) in the title table.");
    clp.RegisterStringOption("filterAdapterOnly", &adapterGffFileName,
            "Use this option to remove reads which can only map to adapters " 
            "specified in the GFF file.");
    clp.SetExamples(
            "Because SAM has optional tags that have different meanings"
            " in different programs, careful usage is required in order "
            "to have proper output.  The \"xs\" tag in bwa-sw is used to "
            "show the suboptimal score, but in PacBio SAM (blasr) it is "
            "defined as the start in the query sequence of the alignment.\n"
            "When \"-smrtTitle\" is specified, the xs tag is ignored, but "
            "when it is not specified, the coordinates given by the xs and "
            "xe tags are used to define the interval of a read that is "
            "aligned.  The CIGAR string is relative to this interval.");

    clp.ParseCommandLine(argc, argv);
    filterCriteria.verbosity = verbosity; 

    // Set random number seed. 
    if (seed == 0) {
        srand(time(NULL));
    } else {
        srand(seed);
    }

    // Set hit policy.
    hitPolicy = setHitPolicy(hitPolicyStr);

    // Set score function.
    scoreFunc = setScoreFunction(scoreFuncStr);

    // Set score cutoff and sign
    if (scoreCutoff != INF_INT)
        filterCriteria.SetScoreCutoff((double)scoreCutoff);
    scoreSign = filterCriteria.SetScoreSign(scoreSignInt);

    string errMsg;
    if (not  filterCriteria.MakeSane(errMsg)) {
        cout << errMsg << endl;
        exit(1);
    }

    if (scoreFunc == USERSCORE and scoreMatrixStr == "") {
        cout << "ERROR. Please specify user-defined score matrix using "
             << "-scoreMatrix." << endl;
        exit(1);
    }

    DistanceMatrixScoreFunction<DNASequence, DNASequence> distScoreFn;
    if (scoreMatrixStr != "") {
        if (scoreFunc != USERSCORE) {
            cout << "ERROR. scoreFunc should be 'userscore' if "
                 << "-scoreMatrix is used." << endl;
            exit(1);
        }
        if (StringToScoreMatrix(scoreMatrixStr, distScoreFn.scoreMatrix) == false) {
            cout << "ERROR. The string " << endl
                << scoreMatrixStr << endl
                << "is not a valid format. It should be a quoted, "
                << "space separated string of "  << endl
                << "integer values.  The matrix: " << endl
                << "    A  C  G  T  N" << endl
                << " A  1  2  3  4  5" << endl
                << " C  6  7  8  9 10" << endl
                << " G 11 12 13 14 15" << endl
                << " T 16 17 18 19 20" << endl
                << " N 21 22 23 24 25" << endl
                << " should be specified as \"1 2 3 4 5 6 7 8 9 "
                << "10 11 12 13 14 15 16 17 18 19 " << endl
                << "20 21 22 23 24 25\"" << endl;
            exit(1);
        }
    }

    // Parse hole number ranges. 
    if (holeNumberStr.size() != 0) {
        if (not holeNumberRanges.setRanges(holeNumberStr)) {
            cout << "Could not parse hole number ranges: "
                 << holeNumberStr << "." << endl;
            exit(1);
        } 
    }

    if (scoreFunc == EDITDIST) {
        // Penalty=1 for each mismatch and indel. 
        distScoreFn.InitializeScoreMatrix(EditDistanceMatrix);
        insScore = delScore = 1;
        scoreSign = filterCriteria.SetScoreSign(NEG);
    } else if (scoreFunc == BLASRSCORE) {
        distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
        insScore = delScore = 5;
        scoreSign = filterCriteria.SetScoreSign(NEG);
    } 
    //
    // If ALIGNERSCORE, score matrix does not matter.
    // If USERSCORE, use user-specified scorematrix, insScore 
    // and delScore.
    //
    distScoreFn.ins = insScore;
    distScoreFn.del = delScore;

    ostream * outFilePtr = &cout;
	ofstream outFileStrm;
	if (outFileName != "") {
		CrucialOpen(outFileName, outFileStrm, std::ios::out);
		outFilePtr = &outFileStrm;
	}
    
    GFFFile adapterGffFile;
    if (adapterGffFileName != "")
        adapterGffFile.ReadAll(adapterGffFileName);
    
    SAMReader<SAMFullReferenceSequence, SAMReadGroup, SAMAlignment> samReader;
    FASTAReader fastaReader;

    //
    // Initialize samReader and fastaReader.
    //
    samReader.Initialize(samFileName);
    fastaReader.Initialize(refFileName);

    //
    // Configure the file log.
    //
    string command;
    CommandLineParser::CommandLineToString(argc, argv, command);
    string log = "Filter sam hits.";
    string program = "samFilter";
    string versionString = VERSION;
    AppendPerforceChangelist(PERFORCE_VERSION_STRING, versionString);

    //
    // Read necessary input.
    //
    vector<FASTASequence> references;
    fastaReader.ReadAllSequences(references);

    // If the SAM file is generated by blasr with -titleTable,
    // then references in the SAM are represented by 
    // their corresponding indices in the title table.
    // In that case, we need to convert reference titles in fasta file
    // to their corresponding indices in the title table, such that
    // references in both SAM and fasta files are represented
    // by title table indices and therefore can match.
    if (titleTableName != "") {
        ConvertTitlesToTitleTableIndices(references, titleTableName);
    }
 
    AlignmentSet<SAMFullReferenceSequence, SAMReadGroup, SAMAlignment> alignmentSet;
    vector<string> allHeaders = samReader.ReadHeader(alignmentSet); 

    // Process SAM Header.
    string commandLineString;
    clp.CommandLineToString(argc, argv, commandLineString);
    allHeaders.push_back("@PG\tID:SAMFILTER\tVN:" + versionString + \
                         "\tCL:" + program + " " + commandLineString);
    for (int i = 0; i < allHeaders.size(); i++) {
        outFileStrm << allHeaders[i] << endl;
    }

    //
    // The order of references in vector<FASTASequence> references and
    // AlignmentSet<, , >alignmentSet.references can be different.
    // Rearrange alignmentSet.references such that they are ordered in
    // exactly the same way as vector<FASTASequence> references.
    //
    alignmentSet.RearrangeReferences(references);

    // Map reference name obtained from SAM file to indices
    map<string, int> refNameToIndex;
    for (int i = 0; i < references.size(); i++) {
        string refName = alignmentSet.references[i].GetSequenceName();
        refNameToIndex[refName] = i;
    }

    //
    // Store the alignments.
    //
    SAMAlignment samAlignment;
    int alignIndex = 0; 

    //
    // For 150K, each chip produces about 300M sequences 
    // (not including quality values and etc.).
    // Let's assume that the sam file and reference data can 
    // fit in the memory. 
    // Need to scale for larger sequal data in the future.
    //
    vector<SAMAlignment> allSAMAlignments;
    while (samReader.GetNextAlignment(samAlignment)) {
        if (samAlignment.rName == "*") {
            continue;
        }

        if (parseSmrtTitle and holeNumberStr.size() != 0) {
            string movieName;
            int thisHoleNumber;
            if (not ParsePBIReadName(samAlignment.qName, 
                                     movieName, 
                                     thisHoleNumber)) {
                cout << "ERROR, could not parse SMRT title: "
                     << samAlignment.qName << "." << endl;
                exit(1);
            }
            if (not holeNumberRanges.contains(UInt(thisHoleNumber))) {
                if (verbosity > 0) 
                    cout << thisHoleNumber << " is not in range." << endl; 
                continue;
            }
        }

        if (samAlignment.cigar.find('P') != string::npos) {
            cout << "WARNING. Could not process SAM record with 'P' in "
                 << "its cigar string." << endl;
            continue;
        }

        vector<AlignmentCandidate<> > convertedAlignments;
        SAMAlignmentsToCandidates(samAlignment, 
                references, refNameToIndex,
                convertedAlignments, parseSmrtTitle, false);
        
        if (convertedAlignments.size() > 1) {
            cout << "WARNING. Ignore multiple segments." << endl;
            continue;
        }

        for (int i = 0; i < 1; i++) {
            AlignmentCandidate<> & alignment = convertedAlignments[i];

            ComputeAlignmentStats(alignment, alignment.qAlignedSeq.seq, 
                                  alignment.tAlignedSeq.seq, distScoreFn);
                                  // scoreMatrix, insScore, delScore);
            if (verbosity > 0)  {
                cout << "Aligner's score = "  << samAlignment.as 
                     << ", computed score = " << alignment.score << endl;
            }

            // Check whether this alignment can only map to adapters in 
            // the adapter GFF file.
            if (adapterGffFileName != "" and 
                CheckAdapterOnly(adapterGffFile, alignment, refNameToIndex)) {
                if (verbosity > 0)
                    cout << alignment.qName << " filter adapter only."
                         << endl;
                continue;
            }

            // Assign score to samAlignment.
            samAlignment.score = alignment.score;
            if (scoreFunc == ALIGNERSCORE) 
                samAlignment.score = samAlignment.as;

            if (not filterCriteria.Satisfy(alignment)) {
                if (verbosity > 0)
                    cout << alignment.qName
                         << " does not satisfy filter criteria." << endl;
                continue;
            }
            allSAMAlignments.push_back( samAlignment ); 

            alignment.FreeSubsequences();
        }
        ++alignIndex;
    }

    // Sort all SAM alignments by qName, score and target position.
    if (!isSorted) {
        sort(allSAMAlignments.begin(), allSAMAlignments.end(), 
             byQNameScoreTStart);
    }

    unsigned int groupBegin = 0;
    unsigned int groupEnd = -1;
    vector<SAMAlignment> filteredSAMAlignments;
    while(groupBegin < allSAMAlignments.size()) {
        // Get the next group of SAM alignments which have the same qName
        // from allSAMAlignments[groupBegin ... groupEnd)
        GetNextSAMAlignmentGroup(allSAMAlignments, groupBegin, groupEnd);
        vector<unsigned int> hitIndices = ApplyHitPolicy(
                hitPolicy, allSAMAlignments, groupBegin, groupEnd);
        for(unsigned int i = 0; i < hitIndices.size(); i++) {
            filteredSAMAlignments.push_back(allSAMAlignments[hitIndices[i]]);
        }
        groupBegin = groupEnd;
    }

    sort(filteredSAMAlignments.begin(), filteredSAMAlignments.end(), 
         byRNameQName);

    for(unsigned int i = 0; i < filteredSAMAlignments.size(); i++) {
        filteredSAMAlignments[i].PrintSAMAlignment(outFileStrm);
    }

	if (outFileName != "") {
		outFileStrm.close();
	}
#ifdef USE_GOOGLE_PROFILER
  ProfilerStop();
#endif
    return 0;
}
