#ifndef FILES_FRAGMENT_CCS_ITERATOR_H_
#define FILES_FRAGMENT_CCS_ITERATOR_H_

#include <vector>
#include "../datastructures/reads/RegionTable.h"
#include "CCSIterator.h"

using namespace std;

class FragmentCCSIterator : public CCSIterator {
public:
  RegionTable *regionTablePtr;
  vector<ReadInterval> subreadIntervals;
  vector<int>          readIntervalDirection;
  virtual void Initialize(CCSSequence *_seqPtr, RegionTable *_regionTablePtr) {
    seqPtr         = _seqPtr;
    regionTablePtr = _regionTablePtr;
    curPass = 0;
    subreadIntervals.clear();
    readIntervalDirection.clear();
    /*
       Since this iterator covers all passes, and not just those
       included in the ccs, the the regions need to be loaded.
    */
    CollectSubreadIntervals(*seqPtr, regionTablePtr, subreadIntervals);
    readIntervalDirection.resize(subreadIntervals.size());
    fill(readIntervalDirection.begin(), readIntervalDirection.end(), 2);

    //
    // Assign the read interval directions based on the pass direction
    // for the pass that has a similar start position.  This allows
    // some wiggle although in practice they coordinates of the pass
    // start base and the template should always match up. 
    //
    
    int i, j;
    bool ccsSubreadFound = false;
    int firstCCSReadPos;
    int firstCCSReadDir;
    int firstInsertThatIsACCSSubread = -1;
    for (i = 0; i < subreadIntervals.size(); i++) {
      for (j = 0; j < seqPtr->passStartBase.size(); j++) {
        if (abs(((int)subreadIntervals[i].start)  - ((int)seqPtr->passStartBase[j])) < 10) {
          readIntervalDirection[i] = seqPtr->passDirection[j];
          break;
        }
      }
    }
    
    int firstAssignedSubread = 0;
    while (firstAssignedSubread < subreadIntervals.size() and 
           readIntervalDirection[firstAssignedSubread] == 2) { 
      firstAssignedSubread++; 
    }
    
    int lastAssignedSubread = subreadIntervals.size();
    while (lastAssignedSubread > 0 and readIntervalDirection[lastAssignedSubread-1] == 2) {
      lastAssignedSubread--;
    }
    
    bool allFullInsertsAreSubreads = true;
    for (i = firstAssignedSubread; i < lastAssignedSubread; i++) {
      if (readIntervalDirection[i] == 2) {
        // Find a read interval whose direction is unknown.
        // This should happen very rarely, just assign 0 as 
        // its direction for now.
        readIntervalDirection[i] = 0;
        allFullInsertsAreSubreads = false;
        break;
      }
    }

    // if (allFullInsertsAreSubreads == false) {
    //    cerr  << "WARNING, there are ccs reads that do not include all full length inserts." << endl;
    // }
    
    if (firstAssignedSubread < subreadIntervals.size() and subreadIntervals.size() > 0) {
      int curSubreadDir = readIntervalDirection[firstAssignedSubread];
            assert(curSubreadDir == 0 or curSubreadDir == 1);
      for (i = firstAssignedSubread - 1; i >= 0; i--) {
        curSubreadDir = (curSubreadDir==0)?1:0;
        readIntervalDirection[i] = curSubreadDir;
      }
    }
    if (lastAssignedSubread > 0 and 
        lastAssignedSubread < subreadIntervals.size() and subreadIntervals.size() > 0){ 
      int curSubreadDir = readIntervalDirection[lastAssignedSubread-1];
      assert(curSubreadDir == 0 or curSubreadDir == 1);
      for (i = lastAssignedSubread; i < subreadIntervals.size(); i++) {
        curSubreadDir = ! curSubreadDir;
        readIntervalDirection[i] = curSubreadDir;
      }
    }
    numPasses = subreadIntervals.size();
  }

  virtual int GetNext(int &direction, int &startBase, int &numBases) {
    if (curPass >= subreadIntervals.size()) {
      return 0;
    }
    direction = int(readIntervalDirection[curPass]);
    startBase = int(subreadIntervals[curPass].start);
    numBases  = int(subreadIntervals[curPass].end - subreadIntervals[curPass].start);
    ++curPass;
    return 1;
  }
};

#endif
