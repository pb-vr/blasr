#ifndef FASTA_SEQUENCE_H_
#define FASTA_SEQUENCE_H_

#include <string>
#include <string.h>
#include "Enumerations.h"
#include "DNASequence.h"
#include "datastructures/reads/ZMWGroupEntry.h"
#include <stdint.h>
using namespace std;
//
// NO proteins for now.
class FASTASequence : public DNASequence {
private:
    bool deleteTitleOnExit;

public:
    char *title;
    int titleLength;
    void PrintSeq(ostream &out, int lineLength = 50, char delim='>') {
        out << delim << title <<endl;
        ((DNASequence*)this)->PrintSeq(out, lineLength);
    }
    int GetStorageSize() {
        if (!title) 
            return DNASequence::GetStorageSize();
        return strlen(title) + DNASequence::GetStorageSize();
    }

    FASTASequence() : DNASequence() {
        title = NULL;
        titleLength = 0;
        deleteTitleOnExit = false; 
        // If deleteTitleOnExist is false, whether to delete title
        // or not should depend on deleteOnExit; otherwise, delete title
        // regardless of deleteOnExit.
    }

    string GetName() {
        string name;
        int i;
        for (i = 0; i < titleLength; i++) {
            if (title[i] != ' ' and
                    title[i] != '\t' and
                    title[i] != '\n' and
                    title[i] != '\r') {
                name.push_back(title[i]);
            }
            else {
                break;
            }
        }
        return name;
    }

	//
	// Define  some no-ops to satisfy instantiating templates that
	// expect these to exist.
	//
	bool StoreHoleNumber(int holeNumber) {return false;}
	bool StoreHoleStatus(unsigned char holeStatus) {return false;}
	bool StorePlatformId(PlatformId platform) {return false;}
	bool StoreZMWData(ZMWGroupEntry &data) { return false;}
	bool GetHoleNumber (int &holeNumberP) {
		//
		// There is no notion of a hole number for a fasta sequence.
		//
		return false;
	}

	bool StoreXY(int16_t xy[]) {return false;};

	bool GetXY(int xyP[]) {
		//
		// Although the xyP is stored in the fasta title for astro reads
		// this class is more general than an astro read, so do not assume 
		// that it may be found in the title.
		//
		// So, this function is effectively a noop.
		//
		xyP[0] = xyP[1] = 0;
		return false;
	}

    void ShallowCopy(const FASTASequence &rhs) {
        CheckBeforeCopyOrReference(rhs, "FASTASequence");
        FASTASequence::Free();

        ((DNASequence*)this)->ShallowCopy(rhs);

        title = rhs.title;
        titleLength = rhs.titleLength;
        deleteTitleOnExit = false;
    }

	string GetTitle() const {
		return string(title);
	}

    void CopyTitle(const char* str, int strlen) {
        FASTASequence::DeleteTitle();

        // No segfault when str is NULL;
        if (str == NULL) {
            title = NULL;
            titleLength = 0;
        } else {
            title = new char[strlen+1];
            memcpy(title, str, strlen);
            titleLength = strlen;
            title[titleLength] = '\0';
        }

        // In some cases, (e.g., when ReferenceSubstring and CopyTitle
        // are called together), this Sequence may only have control over
        // title but not seq.
        deleteTitleOnExit = true;
    }

	void CopyTitle(string str) {
        FASTASequence::CopyTitle(str.c_str(), str.size());
	}

    // Delete title if this FASTASequence is under control or 
    // only title is under control.
    void DeleteTitle() {
        if (deleteOnExit or deleteTitleOnExit) {
            if (title != NULL) {
                delete[] title;
            }
        } // otherwise, title is controlled by another obj
        title = NULL;
        titleLength = 0;
        deleteTitleOnExit = false;
    }

	void GetFASTATitle(string& fastaTitle) {
		// look for the first space, and return the string until there.
		int i;
		for (i = 0; i < titleLength; i++ ){
			if (title[i] == ' ' or
					title[i] == '\t') {
				break;
			}
		}
		fastaTitle.assign(title, i);
	}

    // Copy rhs.seq[readStart:readEnd] to this Sequence.
    void CopySubsequence(FASTASequence &rhs, int readStart, int readEnd) {
        CheckBeforeCopyOrReference(rhs, "FASTASequence");

        // Free before copying anything 
        FASTASequence::Free();

        if (readEnd == -1) {
            readEnd = rhs.length;
        }

        if (readEnd > readStart) {
            length = readEnd - readStart;
            DNASequence::Copy(rhs, readStart, length);
        }
        else {
            seq = NULL;
            length = 0;
            deleteOnExit = true;
        }
        FASTASequence::CopyTitle(rhs.title);
    }

    void AppendToTitle(string str) {
        int newLength = titleLength + str.size() + 1;
        if (newLength == 0) {
            DeleteTitle();
            return;
        }

        char *tmpTitle = new char[newLength];
        memcpy(tmpTitle, title, titleLength);
        memcpy(&tmpTitle[titleLength], str.c_str(), str.size());
        tmpTitle[newLength-1] = '\0';
        delete[] title;
        title = tmpTitle;
        titleLength = newLength;
        deleteTitleOnExit = true;
    }

    void Assign(FASTASequence &rhs) {
        *this = (FASTASequence&)rhs;
    }

    // Create a reverse complement FASTASequence of *this and assign to rhs.
    void MakeRC(FASTASequence &rhs, DNALength rhsPos=0, DNALength rhsLength=0) {
        DNASequence::MakeRC((DNASequence&) rhs, rhsPos, rhsLength);
        if (title != NULL) {
            ((FASTASequence&)rhs).CopyTitle(title);
        }
    }

    void ReverseComplementSelf() {
        DNALength i;
        for (i = 0; i < length/2 + length % 2; i++) {
            char c = seq[i];
            seq[i] = ReverseComplementNuc[seq[length - i - 1]];
            seq[length - i - 1] = ReverseComplementNuc[c];
        }
    }

    void operator=(const FASTASequence &rhs) {
        CheckBeforeCopyOrReference(rhs, "FASTASequence");

        // Free before copying anything 
        FASTASequence::Free();

        // Copy seq from rhs
        ((DNASequence*)this)->Copy((DNASequence&)rhs);

        assert(deleteOnExit);

        // Copy title from rhs
        FASTASequence::CopyTitle(rhs.title, rhs.titleLength);

        assert(deleteOnExit);
    }
	
    void Copy(const FASTASequence &rhs) {
        *this = (FASTASequence&)rhs;
    }

    void Free() {
        // Delete title if title is under control, reset deleteTitleOnExit.
        FASTASequence::DeleteTitle();

        // Delete seq if under control, reset deleteOnExit.
        // Don't call Free() before calling DeleteTitle().
        DNASequence::Free();
    }
};


#endif
