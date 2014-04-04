#ifndef UTILS_FILE_OF_FILE_NAMES_H_
#define UTILS_FILE_OF_FILE_NAMES_H_
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include "../utils.h"
#include "../data/hdf/HDFNewBasReader.h"

class FileOfFileNames {
public:

	static void StoreFileOrFileList(string fileName, vector<string> &fofnList) {
        vector<string> tmpList;
		if (IsFOFN(fileName)) {
			FOFNToList(fileName, tmpList);
		}
		else {
			tmpList.push_back(fileName);
		}
        for (int i = 0; i < int(tmpList.size()); i++) {
            if (FileOfFileNames::IsFOFN(tmpList[i])) {
                cout << "ERROR. Nested File of File Names are not allowed. " 
                     << endl;
                exit(1);
            } else if (FileOfFileNames::IsBasH5(tmpList[i])) {
                vector<string> baxFNs = FileOfFileNames::Bas2Bax(tmpList[i]);
                fofnList.insert(fofnList.end(), baxFNs.begin(), baxFNs.end());
            } else {
                fofnList.push_back(tmpList[i]);
            }
        }
	}

	static void FOFNToList(string &fofnFileName, vector<string> &fofnList) {
		ifstream fofnIn;
		CrucialOpen(fofnFileName, fofnIn);
		while(fofnIn) {
			string name;
			getline(fofnIn, name);
			if (name.size() > 0) {
				fofnList.push_back(name);
			}
		}
	}

	static bool IsFOFN(string &fileName) {
		string::size_type dotPos = fileName.rfind(".");
		if (dotPos != string::npos) {
			string extension;
			extension.assign(fileName, dotPos+1, fileName.size() - (dotPos+1));
			if (extension == "fofn") {
					return true;
				}
		}
		return false;
	}

    static bool IsBasH5(string & fileName) {
        // Return true if file ends with bas.h5
        if (fileName.size() > 6) {
            if (fileName.rfind("bas.h5") == fileName.size() - 6) {
                return true;
            }
        }
        return false;
    }

    static vector<string> Bas2Bax(string & basFN) {
        // There are two types of bas.h5 files. 
        // Before SMRT 2.0, bas.h5 files contain all the /PulseData data, 
        // in this case, return the bas.h5.
        // After SMRT 2.0, bas.h5 files have been changed to only contain 
        // paths to bax.h5 files (in the /MultiPart/Parts group), while
        // all base calls and QVs are in bax.h5 files. In this case, 
        // return path to the bax.h5 files. Assumption is that bax.h5 
        // files are in the same directory as bas.h5 file.
        vector<string> baxFNs; 
        HDFNewBasReader reader;
        if (reader.Initialize(basFN) != 0) {
            baxFNs = reader.GetBaxFileNames();
        } else {
            baxFNs.push_back(basFN);
        }
        reader.Close();
        return baxFNs;
    }

    static int ExpandFileNameList(vector<string> &fileNames) {
        int rfn;
        vector<string> expandedFileNames;
        for (rfn = 0; rfn < fileNames.size(); rfn++) {
            vector<string> tmpList; 
	        FileOfFileNames::StoreFileOrFileList(fileNames[rfn], tmpList);
            expandedFileNames.insert(expandedFileNames.end(), 
                                     tmpList.begin(), tmpList.end());
        }
        fileNames = expandedFileNames;
        return fileNames.size();
    }

};

#endif
