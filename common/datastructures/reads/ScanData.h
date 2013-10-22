#ifndef DATASTRUCTURES_READS_SCAN_DATA_H_
#define DATASTRUCTURES_READS_SCAN_DATA_H_

#include "Enumerations.h"
#include <map>
#include <string>
using namespace std;

class ScanData {
 public:
	PlatformId platformId;
	float frameRate;
	unsigned int numFrames;
	string movieName, runCode;
	string whenStarted;
    map<char, int> baseMap;
	string GetMovieName() {
		return movieName;
	}
    ScanData() {
        platformId = NoPlatform;
        frameRate = numFrames = 0;
        movieName = runCode = whenStarted = "";
        baseMap.clear();
    }
};

#endif
