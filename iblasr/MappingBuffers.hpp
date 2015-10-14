// Author: Mark Chaisson
#ifndef __BLASR_MAPPING_BUFFERS__
#define __BLASR_MAPPING_BUFFERS__

#include <vector>
#include <tuples/DNATuple.hpp>
#include <tuples/TupleList.hpp>
#include <algorithms/alignment/sdp/SDPFragment.hpp>
#include <algorithms/anchoring/BasicEndpoint.hpp>
#include <datastructures/anchoring/ClusterList.hpp>
#include <datastructures/anchoring/MatchPos.hpp>

using namespace std;


//
// Define a list of buffers that are meant to grow to high-water
// marks, and not shrink down past that.   The memory is reused rather
// than having multiple calls to new.
//
class MappingBuffers {
public:
    vector<int> hpInsScoreMat, insScoreMat;
    vector<int> kbandScoreMat;
    vector<Arrow> hpInsPathMat, insPathMat;
    vector<Arrow> kbandPathMat;
    vector<int>   scoreMat;
    vector<Arrow> pathMat;
    vector<int>  affineScoreMat;
    vector<Arrow> affinePathMat;
    vector<ChainedMatchPos> matchPosList;
    vector<ChainedMatchPos> rcMatchPosList;
    vector<BasicEndpoint<ChainedMatchPos> > globalChainEndpointBuffer;
    vector<Fragment> sdpFragmentSet, sdpPrefixFragmentSet, sdpSuffixFragmentSet;
    TupleList<PositionDNATuple> sdpCachedTargetTupleList;
    TupleList<PositionDNATuple> sdpCachedTargetPrefixTupleList;
    TupleList<PositionDNATuple> sdpCachedTargetSuffixTupleList;
    std::vector<int> sdpCachedMaxFragmentChain;
    vector<double> probMat;
    vector<double> optPathProbMat;
    vector<float>  lnSubPValueMat;
    vector<float>  lnInsPValueMat;
    vector<float>  lnDelPValueMat;
    vector<float>  lnMatchPValueMat;
    vector<int>    clusterNumBases;
    ClusterList    clusterList;
    ClusterList    revStrandClusterList;

    void Reset(void);
};


inline void MappingBuffers::Reset(void) {
    vector<int>().swap(hpInsScoreMat);
    vector<int>().swap(insScoreMat);
    vector<int>().swap(kbandScoreMat);
    vector<Arrow>().swap(hpInsPathMat);
    vector<Arrow>().swap(insPathMat);
    vector<Arrow>().swap(kbandPathMat);
    vector<int>().swap(scoreMat);
    vector<Arrow>().swap(pathMat);
    vector<ChainedMatchPos>().swap(matchPosList);
    vector<ChainedMatchPos>().swap(rcMatchPosList);
    vector<BasicEndpoint<ChainedMatchPos> >().swap(globalChainEndpointBuffer);
    vector<Fragment>().swap(sdpFragmentSet);
    vector<Fragment>().swap(sdpPrefixFragmentSet);
    vector<Fragment>().swap(sdpSuffixFragmentSet);
    sdpCachedTargetTupleList.Reset();
    sdpCachedTargetPrefixTupleList.Reset();
    sdpCachedTargetSuffixTupleList.Reset();
    vector<int>().swap(sdpCachedMaxFragmentChain);
    vector<double>().swap(probMat);
    vector<double>().swap(optPathProbMat);
    vector<float>().swap(lnSubPValueMat);
    vector<float>().swap(lnInsPValueMat);
    vector<float>().swap(lnDelPValueMat);
    vector<float>().swap(lnMatchPValueMat);
    vector<int>().swap(clusterNumBases);
}

#endif
