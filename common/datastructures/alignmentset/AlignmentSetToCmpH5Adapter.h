#ifndef ALIGNMENT_SET_TO_CMP_H5_ADAPTER_H_
#define ALIGNMENT_SET_TO_CMP_H5_ADAPTER_H_

#include "utils/SMRTReadUtils.h"
#include "data/hdf/HDFCmpFile.h"
#include "AlignmentSet.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/ByteAlignment.h"
#include "algorithms/alignment/ScoreMatrices.h"

class RefGroupNameId {
 public:
  string name;
  unsigned int id;
  RefGroupNameId() {
    id = 0;
    name = "";
  }

  RefGroupNameId(string n, unsigned int i) {
    name = n; id = i;
  }

  RefGroupNameId(const RefGroupNameId &rhs) {
    name = rhs.name;
    id   = rhs.id;
  }
};

template<typename T_CmpFile>
class AlignmentSetToCmpH5Adapter {
public:
  // Map reference name to reference group (/RefGroup) name and ID.
  map<string, RefGroupNameId> refNameToRefGroupNameandId;
  map<string, unsigned int> knownMovies;
  map<string, unsigned int> knownPaths;
  unsigned int numAlignments;
  map<string, int> refNameToRefInfoIndex;
 
  void Initialize() {
      numAlignments = 0;
  }

  template<typename T_Reference>
  void StoreReferenceInfo(vector<T_Reference> &references, T_CmpFile &cmpFile) {
    for (int r = 0; r < references.size(); r++) {
      string sequenceName, md5;
      sequenceName = references[r].GetSequenceName();
      md5 = references[r].GetMD5();
      unsigned int length = references[r].GetLength();

      // Add this reference to /RefInfo. 
      // Don't create /ref0000x and register it in /RefGroup at this point, 
      // because this reference may not map to any alignments at all. 
      cmpFile.AddRefInfo(sequenceName, length, md5);

      // Update refNameToRefInfoIndex
      if (refNameToRefInfoIndex.find(sequenceName) != refNameToRefInfoIndex.end()) {
        cout << "ERROR. Reference name " << sequenceName 
             << " is not unique." << endl;
        exit(1);
      }
      refNameToRefInfoIndex[sequenceName] = r;
    }
  }

  unsigned int StoreMovieInfo(string movieName, T_CmpFile &cmpFile) {
    map<string, unsigned int>::iterator mapIt;
    mapIt = knownMovies.find(movieName);
    if (mapIt != knownMovies.end()) {
      return mapIt->second;
    }
    else {
      unsigned int id = cmpFile.movieInfoGroup.AddMovie(movieName);
      knownMovies[movieName] = id;
      return id;
    }
  }


  // Given a reference name, find whether there exists a refGroup 
  // (e.g. /ref000001) associated with it. 
  // If not, create a refGroup and update refNameToRefGroupNameandId. 
  // Finally, return its associated refGroup ID.
  unsigned int StoreRefGroup(string refName, T_CmpFile & cmpFile) {
      // Find out whether there is a refGroup associated with refName.
      map<string, RefGroupNameId>::iterator mapIt;
      mapIt = refNameToRefGroupNameandId.find(refName);
      if (mapIt != refNameToRefGroupNameandId.end()) {
        // An existing refGroup is associated with this refName. 
        return mapIt->second.id;
      } else {
        // No refGroup is associated with refName, create one.
        int refInfoIndex = refNameToRefInfoIndex[refName];
        unsigned int refInfoId = refInfoIndex + 1;

        string refGroupName;
        unsigned int refGroupId = cmpFile.AddRefGroup(refName, refInfoId, refGroupName);  

        // Update refNameToRefGroupNameandId.
        refNameToRefGroupNameandId[refName] = RefGroupNameId(refGroupName, refGroupId);
        return refGroupId;
     }
  }

  unsigned int StorePath(string & path, T_CmpFile &cmpFile) {
    if (knownPaths.find(path) != knownPaths.end()) {
      return knownPaths[path];
    }
    else {
      unsigned int id = cmpFile.alnGroupGroup.AddPath(path);
      knownPaths[path] = id;
      return id;
    }
  }

  void RemoveGapsAtEndOfAlignment(AlignmentCandidate<> &alignment) {
    int numEndDel = 0, numEndIns = 0;
    if (alignment.gaps.size() > 0) {
      int lastGap = alignment.gaps.size() - 1;
      int g;
      for (g = 0; g < alignment.gaps[lastGap].size(); g++) {
        if (alignment.gaps[lastGap][g].seq == Gap::Target) {
          numEndIns += alignment.gaps[lastGap][g].length;
        }
        else if (alignment.gaps[lastGap][g].seq == Gap::Query) {
          numEndDel += alignment.gaps[lastGap][g].length;
        }
      }
    }
    alignment.qAlignedSeqLength -= numEndIns;
    alignment.tAlignedSeqLength -= numEndDel;
  }

  void StoreAlignmentCandidate(AlignmentCandidate<> &alignment, 
                               int alnSegment,
                               T_CmpFile &cmpFile, int moleculeNumber = -1) {
    //
    // Find out where the movie is going to get stored.
    //
    string movieName;
    int holeNumber = 0;
    bool nameParsedProperly;
    
    nameParsedProperly = ParsePBIReadName(alignment.qName, movieName, holeNumber);
    if (!nameParsedProperly) {
      cout <<"ERROR. Attempting to store a read with name " 
           << alignment.qName << " that does not " << endl
           << "appear to be a PacBio read." << endl;
      exit(1);
    }
  
    unsigned int movieId = StoreMovieInfo(movieName, cmpFile);

    // Check whether the reference is in /RefInfo.
    map<string, int>::iterator mapIt;
    mapIt = refNameToRefInfoIndex.find(alignment.tName);
    if (mapIt == refNameToRefInfoIndex.end()) {
      cout << "ERROR. The reference name " << alignment.tName 
           << " was not found in the list of references." << endl;
      cout << "Perhaps a different reference file was aligned to than " << endl
           << "what was provided for SAM conversion. " << endl;
      exit(1);
    } 

    // Store refGroup
    unsigned int refGroupId = StoreRefGroup(alignment.tName, cmpFile);
    string refGroupName = refNameToRefGroupNameandId[alignment.tName].name; 
    assert(refGroupId  == refNameToRefGroupNameandId[alignment.tName].id);

    if (cmpFile.refGroupIdToArrayIndex.find(refGroupId) == cmpFile.refGroupIdToArrayIndex.end()) {
      cout << "ERROR. The reference ID is not indexed. " 
           << "This is an internal inconsistency." << endl;
      exit(1);
    }

    int    refGroupIndex= cmpFile.refGroupIdToArrayIndex[refGroupId];
    assert(refGroupIndex + 1 == refGroupId);

    string path = "/" + refGroupName + "/" + movieName;
    unsigned int pathId = StorePath(path, cmpFile);
    int pathIndex = pathId - 1;

    vector<unsigned int> alnIndex;
    alnIndex.resize(22);

    RemoveGapsAtEndOfAlignment(alignment);
  
    /*
     * Store the alignment string
     */
    vector<unsigned char> byteAlignment;
    AlignmentToByteAlignment(alignment, 
                             alignment.qAlignedSeq, alignment.tAlignedSeq,
                             byteAlignment);

    unsigned int offsetBegin, offsetEnd;
    cmpFile.StoreAlnArray(byteAlignment, alignment.tName, movieName, offsetBegin, offsetEnd);

    numAlignments++;
    /*    EditDistanceMatrix scoreMat;
    */
    int tmpMatrix[5][5];
    // the 5,5 are indel penalties that do not matter since the score is not stored.
    ComputeAlignmentStats(alignment, alignment.qAlignedSeq.seq, alignment.tAlignedSeq.seq, tmpMatrix, 5, 5);

    /*
      The current AlnIndex column names:
      (0): "AlnID", "AlnGroupID", "MovieID", "RefGroupID", "tStart",
      (5): "tEnd", "RCRefStrand", "HoleNumber", "SetNumber",
      (9): "StrobeNumber", "MoleculeID", "rStart", "rEnd", "MapQV", "nM",
      (15): "nMM", "nIns", "nDel", "Offset_begin", "Offset_end",
      (20): "nBackRead", "nReadOverlap"
    */
    if (moleculeNumber == -1) {
      moleculeNumber = holeNumber * movieId;
    }
    alnIndex[0]  = numAlignments;  // AlnId
    alnIndex[1]  = pathId;        // AlnGroupID
    alnIndex[2]  = movieId;    // MovieID
    alnIndex[3]  = refGroupId; // RefGroupID
    alnIndex[4]  = alignment.tAlignedSeqPos; // tStart
    alnIndex[5]  = alignment.tAlignedSeqPos +  alignment.tAlignedSeqLength; // tEnd
    alnIndex[6]  = alignment.tStrand; // RCRefStrand
    alnIndex[7]  = holeNumber;
    alnIndex[8]  = 0; // SET NUMBER -- parse later!!!!
    alnIndex[9]  = alnSegment; // strobenumber
    alnIndex[10] = moleculeNumber;
    alnIndex[11] = alignment.qAlignedSeqPos; 
    alnIndex[12] = alignment.qAlignedSeqPos + alignment.qAlignedSeqLength;
    alnIndex[13] = alignment.mapQV;
    alnIndex[14] = alignment.nMatch;
    alnIndex[15] = alignment.nMismatch;
    alnIndex[16] = alignment.nIns;
    alnIndex[17] = alignment.nDel;
    alnIndex[18] = offsetBegin;
    alnIndex[19] = offsetEnd;
    alnIndex[20] = 0;
    alnIndex[21] = 0;
    cmpFile.alnInfoGroup.WriteAlnIndex(alnIndex);
  }

  void StoreAlignmentCandidateList(vector<AlignmentCandidate<> > &alignments, T_CmpFile &cmpFile,  int moleculeNumber=-1) {
    int a;
    for (a = 0; a < alignments.size(); a++) {
      StoreAlignmentCandidate(alignments[a], a, cmpFile, moleculeNumber);
    }
  }

  void StoreAlignmentCandidate(AlignmentCandidate<> alignment, 
                               T_CmpFile &cmpFile) {
    StoreAlignmentCandidate(alignment, 0, cmpFile);
  }

};



#endif
