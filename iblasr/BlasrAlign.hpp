// Author: Mark Chaisson
#pragma once

#include "BlasrHeaders.h"
#include "BlasrMiscs.hpp"

//------------------MAP READS---------------------------------//
template<typename T_Sequence, typename T_RefSequence, typename T_SuffixArray, typename T_TupleCountTable>
void MapRead(T_Sequence &read, T_Sequence &readRC, T_RefSequence &genome,
             T_SuffixArray &sarray,
             BWT &bwt,
             SeqBoundaryFtr<FASTQSequence> &seqBoundary,
             T_TupleCountTable &ct,
             SequenceIndexDatabase<FASTQSequence> &seqdb,
             MappingParameters &params,
             MappingMetrics    &metrics,
             vector<T_AlignmentCandidate*> &alignmentPtrs,
             MappingBuffers &mappingBuffers,
             MappingIPC *mapData,
             MappingSemaphores & semaphores);

template<typename T_Sequence>
void MapRead(T_Sequence &read, T_Sequence &readRC,
             vector<T_AlignmentCandidate*> &alignmentPtrs,
             MappingBuffers &mappingBuffers,
             MappingIPC *mapData,
             MappingSemaphores & semaphores);

/*
void MapReads(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData);
*/

//------------------MAKE ALIGNMENTS---------------------------//
template<typename T_TargetSequence, typename T_QuerySequence, typename TDBSequence>
void AlignIntervals(T_TargetSequence &genome, T_QuerySequence &read, T_QuerySequence &rcRead,
                    WeightedIntervalSet &weightedIntervals,
                    int mutationCostMatrix[][5],
                    int ins, int del, int sdpTupleSize,
                    int useSeqDB, SequenceIndexDatabase<TDBSequence> &seqDB,
                    vector<T_AlignmentCandidate*> &alignments,
                    MappingParameters &params,
                    MappingBuffers &mappingBuffers,
                    int procId=0);

template<typename T_RefSequence, typename T_Sequence>
void PairwiseLocalAlign(T_Sequence &qSeq, T_RefSequence &tSeq,
                        int k,
                        MappingParameters &params, T_AlignmentCandidate &alignment,
                        MappingBuffers &mappingBuffers,
                        AlignmentType alignType=Global);

// Extend target aligned sequence of the input alignement to both ends
// by flankSize bases. Update alignment->tAlignedSeqPos,
// alignment->tAlignedSeqLength and alignment->tAlignedSeq.
void FlankTAlignedSeq(T_AlignmentCandidate * alignment,
                      SequenceIndexDatabase<FASTQSequence> &seqdb,
                      DNASequence & genome,
                      int flankSize);

// Align a subread of a SMRT sequence to target sequence of an alignment.
// Input:
//   subread         - a subread of a SMRT sequence.
//   unrolledRead    - the full SMRT sequence.
//   alignment       - an alignment.
//   passDirection   - whether or not the subread has the
//                     same direction as query of the alignment.
//                     0 = true, 1 = false.
//   subreadInterval - [start, end) interval of the subread in the
//                     SMRT read.
//   subreadIndex    - index of the subread in allReadAlignments.
//   params          - mapping paramters.
// Output:
//   allReadAlignments - where the sequence and alignments of the
//                       subread are saved.
//   threadOut         - an out stream for debugging the current thread.
void AlignSubreadToAlignmentTarget(ReadAlignments & allReadAlignments,
        SMRTSequence & subread, SMRTSequence & unrolledRead,
        T_AlignmentCandidate * alignment,
        int passDirection, ReadInterval & subreadInterval,
        int subreadIndex,
        MappingParameters & params,
        MappingBuffers & mappingBuffers,
        ostream & threadOut);

#include "BlasrAlignImpl.hpp"
