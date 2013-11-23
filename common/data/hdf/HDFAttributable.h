#ifndef DATA_HDF_HDF_ATTRIBUTABLE_H_
#define DATA_HDF_HDF_ATTRIBUTABLE_H_

#include "H5Cpp.h"
#include <string>
#include <vector>
#include <assert.h>

using namespace std;
using namespace H5;

#ifdef __H5Location_H
void CallStoreAttributeName(H5Location &obj, string attrName, void *attrListPtr);
#else
void CallStoreAttributeName(H5Object &obj, string attrName, void *attrListPtr);
#endif

class HDFAttributable {
 public:
	vector<string> attributeNameList;

 private:

#ifdef __H5Location_H
	void StoreAttributeNames(H5Location &thisobject, vector<string> &attributeNames) {
#else
	void StoreAttributeNames(H5Object &thisobject, vector<string> &attributeNames) {
#endif
		void *destAndData[2];
		int nAttr = thisobject.getNumAttrs();
		unsigned int bounds[2];
		bounds[0] = 0;
		bounds[1] = nAttr;
		attributeNameList.clear();
		thisobject.iterateAttrs(&CallStoreAttributeName, 
								bounds, (void*) &attributeNames);
	}

 public:
  virtual H5Object* GetObject() {
    return NULL;
  }

	int ContainsAttribute(string attributeName) {
		int i;
    vector<string> tmpAttributeNames;
#ifdef __H5Location_H
    H5Location *obj = GetObject();
#else
    H5Object *obj = GetObject();
#endif

    assert(obj != NULL);
    StoreAttributeNames(*obj, tmpAttributeNames);
		for (i = 0; i < tmpAttributeNames.size(); i++) {
			if (tmpAttributeNames[i] == attributeName) return true;
		}
		return false;
	}

};

#ifdef __H5Location_H
void CallStoreAttributeName(H5Location &obj, string attrName, void *attrList){ 
#else
void CallStoreAttributeName(H5Object &obj, string attrName, void *attrList){ 
#endif
	((vector<string>*)attrList)->push_back(attrName);
 }

#endif
