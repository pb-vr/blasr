/*
 * ============================================================================
 *
 *       Filename:  HDFUtils.h
 *
 *    Description:  Implement HDF utils.
 *
 *        Version:  1.0
 *        Created:  11/04/2013 10:13:11 AM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Yuan Li (yli), yli@pacificbiosciences.com
 *        Company:  Pacific Biosciences
 *
 * ============================================================================
 */
#include "HDFFile.h"
#include "HDFScanDataReader.h"
#include "HDFRegionTableReader.h"
#include "datastructures/reads/RegionTable.h"
#include <vector>
#include <string>

using namespace std;

// Given a PacBio (pls/plx/bas/bax/ccs/rgn).h5 file, which contains its movie 
// name in group /ScanData/RunInfo attribute MovieName, return its' movie name
string GetH5MovieName(string fileName) {
    HDFScanDataReader reader;
    return reader.GetMovieName_and_Close(fileName);
}

// Given a vector of h5 files, return their movie names.
vector<string> GetH5MovieNames(const vector<string> & fileNames) {
    vector<string> ret;
    for (int i = 0 ; i < fileNames.size(); i++) {
        ret.push_back(GetH5MovieName(fileNames[i]));
    }
    return ret;
}

// Given a PacBio rgn.h5 file, return the smallest and largest holeNumber in
// group /PulseData/Regions.
pair<UInt, UInt> GetMinMaxHoleNumber(string fileName, 
    bool isRGN=false) {
    UInt minHole, maxHole;

    if (isRGN) { // is region table
        HDFRegionTableReader rgnReader;
        rgnReader.Initialize(fileName);
        rgnReader.GetMinMaxHoleNumber(minHole, maxHole);
        rgnReader.Close();
    } else { // is bas/bax/pls/plx/ccs.h5
        HDFBasReader basReader;
        basReader.Initialize(fileName);
        vector<UInt> holes;
        basReader.GetMinMaxHoleNumber(minHole, maxHole);
        basReader.Close();
    }
    return make_pair(minHole, maxHole);
}

vector< pair<UInt, UInt> > GetMinMaxHoleNumbers(
    const vector<string> & fileNames, bool isRGN=false) {
    vector< pair<UInt, UInt> > ret;
    for (int i = 0 ; i < fileNames.size(); i++) {
        ret.push_back(GetMinMaxHoleNumber(fileNames[i], isRGN));
    }
    return ret;
}

//
// Pulse files in input.fofn and regions tables in rgn.fofn may not
// match, return mapping from plsFNs indices to rgnFNs indices.
//
// Input : plsFNs - pulse file names in input.fofn, e.g.,
//                  P=(p_0, ..., p_{n-1})
//         rgnFNs - region table file names in rgn.fofn, e.g.,
//                  R=(r_0, ..., p_{n-1})
// Output: mapping from plsFNs indices to rgnFNs indices, e.g.,
//                  M=(m_0, ..., m_{n-1})
//         so that for all i from 0 to n-1,
//                  r_{m_{i}} matches p_i
//
vector<int> MapPls2Rgn(const vector<string> & plsFNs,
        const vector<string> & rgnFNs) {
    if (plsFNs.size() != rgnFNs.size() && rgnFNs.size() != 0) {
        cout << "ERROR, the number of plx/bax.h5 files and the number of "
            << "region tables are not the same." << endl;
        exit(1);
    }

    // Movie names of pulse files in P.
    vector<string> plsMovies = GetH5MovieNames(plsFNs);
    // Movie names of region tables in R.
    vector<string> rgnMovies = GetH5MovieNames(rgnFNs);

    // The first and last hole numbers of pulse files in P.
    vector< pair<UInt, UInt> > plsHoles = GetMinMaxHoleNumbers(plsFNs);
    // The first and last hole numbers of region tables in R.
    vector< pair<UInt, UInt> > rgnHoles = GetMinMaxHoleNumbers(rgnFNs, true);

    vector<int> ret;
    for (int i = 0; i < plsFNs.size(); i++) {
        int j = 0;
        for (; j < rgnFNs.size(); j++) {
            if (plsMovies[i] == rgnMovies[j] and
                plsHoles[i].first <= rgnHoles[j].first and
                plsHoles[i].second >= rgnHoles[j].second) {
                break;
            }
        }
        if (j >= rgnFNs.size()) {
            cout << "ERROR, could not find any region table for file "
                 << plsFNs[i] << " [" << plsHoles[i].first << ", " << plsHoles[i].second 
                 <<"." << endl;
            exit(1);
        }
        ret.push_back(j);
    }
    return ret;
}

