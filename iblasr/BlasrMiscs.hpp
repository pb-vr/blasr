// Author: Mark Chaisson
#ifndef _BLASR_MISCS_HPP_
#define _BLASR_MISCS_HPP_

#include "BlasrHeaders.h"

//-------------------------Fetch Reads----------------------------//
template<typename T_Sequence>
bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader,
                                 MappingParameters &params,
                                 T_Sequence &read,
                                 string & readGroupId,
                                 int & associatedRandInt,
                                 MappingSemaphores & semaphores);

//---------------------MAKE & CHECK READS-------------------------//
//FIXME: move to SMRTSequence
bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence);

//FIXME: Move to SMRTSequence
// Given a SMRT sequence and a subread interval, make the subread.
// Input:
//   smrtRead         - a SMRT sequence
//   subreadInterval  - a subread interval
//   params           - mapping parameters
// Output:
//   subreadSequence - the constructed subread
void MakeSubreadOfInterval(SMRTSequence & subreadSequence,
                           SMRTSequence & smrtRead,
                           ReadInterval & subreadInterval,
                           MappingParameters & params);

//FIXME: Move to SMRTSequence
// Given a SMRT sequence and one of its subreads, make the
// reverse complement of the subread in the coordinate of the
// reverse complement sequence of the SMRT sequence.
// Input:
//   smrtRead          - a SMRT read
//   subreadSequence   - a subread of smrtRead
// Output:
//   subreadSequenceRC - the reverse complement of the subread
//                       in the coordinate of the reverse
//                       complement of the SMRT read.
void MakeSubreadRC(SMRTSequence & subreadSequenceRC,
                   SMRTSequence & subreadSequence,
                   SMRTSequence & smrtRead);

// Make a virtual SMRTSequence (polymerase reads) given all subreads.
// NO QVs will be copied at this point.
void MakeVirtualRead(SMRTSequence & smrtRead,
                     const vector<SMRTSequence> & subreads);

// Construct subreads invervals from subreads
void MakeSubreadIntervals(vector<SMRTSequence> & subreads,
                          vector<ReadInterval> & subreadIntervals);

// Get index of median length interval
int GetIndexOfMedian(const vector<ReadInterval> & subreadIntervals);

//-------------------------MISC-----------------------------------//
int CountZero(unsigned char *ptr, int length);

#include "BlasrMiscsImpl.hpp"

#endif
