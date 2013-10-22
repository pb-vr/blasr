#ifndef DATA_HDF_HDF_SCAN_DATA_WRITER_H_
#define DATA_HDF_HDF_SCAN_DATA_WRITER_H_

#include "HDFFile.h"
#include "HDFGroup.h"
#include "HDFAtom.h"
#include "Enumerations.h"
#include "datastructures/reads/ScanData.h"

class HDFScanDataWriter {
private:
    HDFGroup * rootGroupPtr;
	HDFGroup scanDataGroup;
	HDFGroup acqParamsGroup;
    HDFGroup dyeSetGroup;
	HDFGroup runInfoGroup;

	HDFAtom<string> whenStartedAtom;
	HDFAtom<float> frameRateAtom;
	HDFAtom<unsigned int> numFramesAtom;

    HDFAtom<string> baseMapAtom;
    HDFAtom<unsigned int> numAnalogAtom;

	HDFAtom<string> movieNameAtom;
	HDFAtom<string> runCodeAtom;

	HDFAtom<unsigned int> platformIdAtom;
	HDFAtom<string> platformNameAtom;

    void CreateAcqParamsGroup() {
        if (acqParamsGroup.Initialize(scanDataGroup, "AcqParams") == 0) {
            cout << "ERROR could not create /ScanData/AcqParams." << endl;
            exit(1);
        }
        frameRateAtom.Create(acqParamsGroup.group, "FrameRate");
        numFramesAtom.Create(acqParamsGroup.group, "NumFrames");
        whenStartedAtom.Create(acqParamsGroup.group, "WhenStarted");
    }

    void CreateDyeSetGroup(){
        if (dyeSetGroup.Initialize(scanDataGroup, "DyeSet") == 0) {
            cout << "ERROR could not create /ScanData/DyeSet." << endl;
            exit(1);
        }
        baseMapAtom.Create(dyeSetGroup.group, "BaseMap");
        numAnalogAtom.Create(dyeSetGroup.group, "NumAnalog");
    }

    void CreateRunInfoGroup(){
        if (runInfoGroup.Initialize(scanDataGroup, "RunInfo") == 0) {
            cout << "ERROR, could not create /ScanDta/RunInfo." << endl;
            exit(1);
        }
        movieNameAtom.Create(runInfoGroup.group, "MovieName");
        platformIdAtom.Create(runInfoGroup.group, "PlatformId");
        platformNameAtom.Create(runInfoGroup.group, "PlatformName");
        runCodeAtom.Create(runInfoGroup.group, "RunCode");
    }

public:
	HDFScanDataWriter(HDFFile & _outFile) {
        Initialize(_outFile.rootGroup);
	}
    HDFScanDataWriter(HDFGroup & _rootGroup) {
        Initialize(_rootGroup);
    }
    ~HDFScanDataWriter() { 
		// Assume that closing the hdf file must be done
		// manually and not in a destructor.
	}

    int Initialize(HDFGroup & _rootGroup) {
        rootGroupPtr = &(_rootGroup);
        rootGroupPtr->AddGroup("ScanData"); 
        if (scanDataGroup.Initialize(*(rootGroupPtr), "ScanData") == 0) {
            cout << "ERROR, could not create /ScanData group." << endl;
            exit(1);
        }
        scanDataGroup.AddGroup("AcqParams");
        scanDataGroup.AddGroup("DyeSet");
        scanDataGroup.AddGroup("RunInfo");

        CreateAcqParamsGroup();
        CreateDyeSetGroup();
        CreateRunInfoGroup();
    }

    void Write(ScanData & scanData) {
        WriteFrameRate((scanData.frameRate==0)?
                       (75):(scanData.frameRate));
		WriteNumFrames((scanData.numFrames==0)?
                       (1000000):(scanData.numFrames));
		WriteWhenStarted((scanData.whenStarted.empty())?
                        ("2013-01-01T01:01:01"):(scanData.whenStarted));
        string baseMapStr = BaseMapToStr(scanData.baseMap);
        WriteBaseMap((baseMapStr == "")?("TGAC"):baseMapStr);
        WriteNumAnalog(4);

		WriteMovieName((scanData.movieName.empty()?
                       ("simulated_movie"):scanData.movieName));
    	WriteRunCode((scanData.runCode.empty())?
                      "simulated_runcode":(scanData.runCode));
        WritePlatformId((scanData.platformId==NoPlatform)?
                        (Springfield):(scanData.platformId));
	}

	void WriteFrameRate(float frameRate) {
        // Write /ScanData/AcqParams/FrameRate attribute.
        frameRateAtom.Write(frameRate);
    }

    void WriteNumFrames(unsigned int numFrames) {
        // Write /ScanData/AcqParams/NumFrames attribute.
        numFramesAtom.Write(numFrames);
    }

    void WriteWhenStarted(const string whenStarted) {
        // Write /ScanData/AcqParams/WhenStarted attribute.
		whenStartedAtom.Write(whenStarted);
	}

    string BaseMapToStr(map<char, int> & baseMap) {
        string baseMapStr = ""; //4 dye channels.
        if (not baseMap.empty()) {
            baseMapStr = "    ";
            map<char, int>::iterator it;
            for (it = baseMap.begin(); it != baseMap.end(); ++it){
                if (it->second > 4 or it->second < 0) {
                    cout << "ERROR, there are more than four dye channels."
                        << endl;
                    exit(1);
                }
                baseMapStr[it->second]= it->first;
            }
        }
        return baseMapStr;
    }

    void WriteBaseMap(const string baseMapStr) {
        //Write /ScanData/DyeSet/BaseMap attribute.
        baseMapAtom.Write(baseMapStr);
    }

    void WriteNumAnalog(const unsigned int numAnalog) {
        //Write /ScanData/DyeSet/NumAnalog attribute.
        numAnalogAtom.Write(numAnalog);
    }

    void WritePlatformId(const PlatformId id) {
        //Write /ScanData/RunInfo/Flatform attribute.
        platformIdAtom.Write(id);
        string name = (id == Springfield)?"Springfield":"Astro";
        platformNameAtom.Write(name);
    }

    void WriteMovieName(const string movieName) {
        //Write /ScanData/RunInfo/MovieName attribute.
		movieNameAtom.Write(movieName);
	}

    void WriteRunCode(const string runCode) {
        //Write /ScanData/RunInfo/MovieName attribute.
		runCodeAtom.Write(runCode);
	}

	void Close() {
        // Close /ScanData/AcqParams attributes.
        whenStartedAtom.dataspace.close();
	    frameRateAtom.dataspace.close();
	    numFramesAtom.dataspace.close();

        // Close /ScanData/DyeSet attributes.
        baseMapAtom.dataspace.close();

        // Close /ScanData/RunInfo attributes.
		movieNameAtom.dataspace.close();
		runCodeAtom.dataspace.close();
	    platformIdAtom.dataspace.close();
	    platformNameAtom.dataspace.close();

        // Close /ScanData/AcqParams|DyeSet|RunInfo.
		acqParamsGroup.Close();
        dyeSetGroup.Close();
		runInfoGroup.Close();

        // Close /ScanData
		scanDataGroup.Close();
    }
};

#endif
