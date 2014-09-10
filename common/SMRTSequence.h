#ifndef  SMRT_SEQUENCE_H_
#define  SMRT_SEQUENCE_H_
#include <assert.h>
#include <iostream>
#include <fstream>
#include <string>

#include "NucConversion.h"
#include "FASTQSequence.h"
#include "Enumerations.h"
#include "datastructures/reads/RegionTable.h"
#include "datastructures/reads/ZMWGroupEntry.h"
using namespace std;
#include <iostream>
#include <sstream>
#include "Types.h"

typedef unsigned char Nucleotide;

class SMRTSequence : public FASTQSequence {
public:
  int16_t xy[2];
  int holeNumber;
  ZMWGroupEntry zmwData;
  PlatformId platform;
  HalfWord *preBaseFrames;
  HalfWord *widthInFrames;
  //
  // The following are fields that are read in from the pulse file.
  // Because they are not standard in bas.h5 files, these fields
  // should not be preallocated when resizing a SMRTSequence, and
  // memory should be managed separately.  For now, these fields all
  // have the same length as the number of bases, but this could
  // change so that all pulse values are stored in a SMRTSequence.
  //
  HalfWord *meanSignal, *maxSignal, *midSignal;
  float *classifierQV;
  unsigned int *startFrame;
  int *pulseIndex;
  DNALength lowQualityPrefix, lowQualitySuffix;

  void SetNull() {
    pulseIndex    = NULL;
    preBaseFrames = NULL;
    widthInFrames = NULL;
    xy[0] = 0; xy[1] = 0;
    // These are not allocted by default.
    meanSignal = maxSignal = midSignal = NULL;
    classifierQV = NULL;
    startFrame   = NULL;
    platform     = NoPlatform;
    // By default, allow the entire read.
    lowQualityPrefix = lowQualitySuffix = 0;
  }
  
  SMRTSequence() : FASTQSequence() {
    holeNumber = -1;
    SetNull();
  }

  void Allocate(DNALength length) {
      // Assert *this has no allocated space.
      if (not (seq == NULL && preBaseFrames == NULL &&
               widthInFrames == NULL and pulseIndex == NULL)) {
          cout << "ERROR, trying to double-allocate memory for a SMRTSequence." << endl;
          exit(1);
      }

      FASTQSequence::AllocateRichQualityValues(length);
      seq           = new Nucleotide[length];
      qual.Allocate(length);
      preBaseFrames = new HalfWord[length];
      widthInFrames = new HalfWord[length];
      pulseIndex    = new int[length];
      subreadEnd    = length;
      deleteOnExit  = true;
  }

  void SetSubreadTitle(SMRTSequence &subread, DNALength subreadStart, DNALength  subreadEnd) {
    stringstream titleStream;
    titleStream << title << "/"<< subreadStart << "_" << subreadEnd;
    subread.CopyTitle(titleStream.str());
  }    

  void SetSubreadBoundaries(SMRTSequence &subread, DNALength &subreadStart, int &subreadEnd) {
    if (subreadEnd == -1) {
      subreadEnd = length;
    }
    assert(subreadEnd - subreadStart <= length);
    subread.subreadStart= subreadStart;
    subread.subreadEnd  = subreadEnd;
    SetSubreadTitle(subread, subreadStart, subreadEnd);
  }

  void MakeSubreadAsMasked(SMRTSequence &subread, DNALength subreadStart = 0, int subreadEnd = -1) {
      subread.Free();
      //
      // This creates the entire subread, but masks out the portions
      // that do not correspond to this insert.
      //
      ((SMRTSequence&)subread).Copy(*this);
      SetSubreadBoundaries(subread, subreadStart, subreadEnd);
      DNALength pos;
      for (pos = 0; pos < subreadStart; pos++) { subread.seq[pos] = 'N'; }
      for (pos = subreadEnd; pos < length; pos++) { subread.seq[pos] = 'N'; }
      // This is newly allocated memory, free it on exit.
      assert(subread.deleteOnExit);
  }

  void MakeSubreadAsReference(SMRTSequence &subread, DNALength subreadStart = 0, int subreadEnd = -1) {
      subread.Free();
      //
      // Just create a reference to a substring of this read.  
      //
      ((FASTQSequence)subread).ReferenceSubstring(*this, subreadStart, subreadEnd - subreadStart);
      SetSubreadBoundaries(subread, subreadStart, subreadEnd);
      // The subread references this read, protect the memory.
      assert(not subread.deleteOnExit);
  }

  void Copy(const SMRTSequence &rhs) {
      SMRTSequence::Copy(rhs, 0, rhs.length);
  }

  void Copy(const SMRTSequence &rhs, int rhsPos, int rhsLength) {
      // Sanity check
      CheckBeforeCopyOrReference(rhs, "SMRTSequence");

      // Free this SMRTSequence before copying anything from rhs.
      SMRTSequence::Free();

      FASTQSequence subseq; 
      // subseq.seq is referenced, while seq.title is not, we need to call 
      // subseq.Free() to prevent memory leak.
      ((FASTQSequence&)subseq).ReferenceSubstring((FASTQSequence&)rhs, rhsPos, rhsLength);
      ((FASTQSequence&)subseq).CopyTitle(rhs.title, rhs.titleLength); 

      if (rhs.length == 0) {
          ((FASTQSequence*)this)->Copy(subseq);
          //
          // Make sure that no values of length 0 are allocated by returning here.
          //
      }
      else {
          assert(rhs.seq != seq);
          assert(rhsLength <= rhs.length);
          assert(rhsPos < rhs.length);

          // Copy seq, title and FASTQ QVs from subseq
          ((FASTQSequence*)this)->Copy(subseq); 

          // Copy SMRT QVs
          if (rhs.preBaseFrames != NULL) {
              preBaseFrames = new HalfWord[length];
              memcpy(preBaseFrames, rhs.preBaseFrames, length*sizeof(HalfWord));
          }
          if (rhs.widthInFrames != NULL) {
              widthInFrames = new HalfWord[length];
              memcpy(widthInFrames, rhs.widthInFrames, length*sizeof(HalfWord));
          }
          if (rhs.pulseIndex != NULL) {
              pulseIndex = new int[length];
              memcpy(pulseIndex, rhs.pulseIndex, sizeof(int) * length);
          }
      }

      // Copy other member variables from rhs
      subreadStart = rhs.subreadStart;
      subreadEnd   = rhs.subreadEnd;
      lowQualityPrefix = rhs.lowQualityPrefix;
      lowQualitySuffix = rhs.lowQualitySuffix;
      zmwData = rhs.zmwData;

      assert(deleteOnExit); // should have control over seq and all QVs

      subseq.Free();
  }

  void Print(ostream &out) {
    out << "SMRTSequence for zmw " << zmwData.holeNumber
        << ", [" << subreadStart << ", " << subreadEnd << ")" << endl;
    DNASequence::Print(out);
  }

  SMRTSequence& operator=(const SMRTSequence &rhs) {
      SMRTSequence::Copy(rhs);
      return *this;
  }

  void Free() {
      if (deleteOnExit == true) {
          if (preBaseFrames)  {
              delete[] preBaseFrames;
          }
          if (widthInFrames) {
              delete[] widthInFrames;
          }
          if (pulseIndex) {
              delete[] pulseIndex;
          }
          if (startFrame) {
              delete[] startFrame;
          }
          // meanSignal, maxSignal, midSignal and classifierQV
          // need to be handled separatedly.
      }

      // Reset SMRT QV pointers anyway
      preBaseFrames = NULL;
      widthInFrames = NULL;
      pulseIndex = NULL;
      startFrame = NULL;

      // Reset member variables
      xy[0] = 0; xy[1] = 0;
      lowQualityPrefix = lowQualitySuffix = 0;
      holeNumber = -1;

      // Free seq, title and FASTQ QVs, also reset deleteOnExit.
      // Don't call FASTQSequence::Free() before freeing SMRT QVs.
      FASTQSequence::Free();
  }

  bool StoreXY(int16_t xyP[]) {
    xy[0] = xyP[0];
    xy[1] = xyP[1];
    return true;
  }

  bool StorePlatformId(PlatformId pid) {
      platform = pid;
      return true;
  }

  bool StoreHoleNumber(int holeNumberP){ 
    zmwData.holeNumber = holeNumber = holeNumberP;
    return true;
  }
  
  bool StoreHoleStatus(unsigned int s) {
    zmwData.holeStatus = s;
    return true;
  }
  
  bool StoreZMWData(ZMWGroupEntry &data) {
    zmwData = data;
    return true;
  }

  bool GetXY(int xyP[]) {
    xyP[0] = xy[0];
    xyP[1] = xy[1];
    return true;
  }

  bool GetHoleNumber(int& holeNumberP) {
    holeNumberP = holeNumber;
    return true;
  }
};

#endif
