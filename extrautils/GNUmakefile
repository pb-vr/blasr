.PHONY=all cramtests

ROOT:=..
include ${ROOT}/q.mk

ALL_CXXFLAGS := \
		   ${CXXFLAGS}

EXE = sa2bwt bwt2sa alchemy excrep evolve bsdb simpleShredder swMatcher \
  samodify sals printTupleCountTable cmpH5StoreQualityByContext ccsh5tobam

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH


all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${ALL_CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

sa2bwt: SuffixArrayToBWT.cpp
bwt2sa: BwtToSuffixArray.cpp
alchemy: BasH5Simulator.cpp
excrep: ExciseRepeats.cpp
evolve: Evolve.cpp
bsdb: BuildSequenceDB.cpp
simpleShredder: SimpleShredder.cpp
swMatcher: SWMatcher.cpp
samodify: SAModify.cpp
sals: SALS.cpp
printTupleCountTable: PrintTupleCountTable.cpp
cmpH5StoreQualityByContext: StoreQualityByContextFromCmpH5.cpp
ccsh5tobam: CCSH5ToBam.cpp

CTESTS := $(wildcard ctest/*.t)

cramtests: ${EXE}
	cram -v --shell=/bin/bash ${CTESTS}

clean: 
	@rm -f ${EXE}
	@rm -f *.d *.o
