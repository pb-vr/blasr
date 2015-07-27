.PHONY=all cramtests

THISDIR:=$(dir $(lastword $(MAKEFILE_LIST)))
ROOT:=${THISDIR}/..
include ${CURDIR}/../defines.mk
include ${ROOT}/rules.mk

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}

EXE = sa2bwt bwt2sa alchemy excrep evolve bsdb simpleShredder swMatcher \
  samodify sals printTupleCountTable cmpH5StoreQualityByContext ccsh5tobam

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH


all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

sa2bwt: SuffixArrayToBWT.o
bwt2sa: BwtToSuffixArray.o
alchemy: BasH5Simulator.o
excrep: ExciseRepeats.o
evolve: Evolve.o
bsdb: BuildSequenceDB.o
simpleShredder: SimpleShredder.o
swMatcher: SWMatcher.o
samodify: SAModify.o
sals: SALS.o
printTupleCountTable: PrintTupleCountTable.o
cmpH5StoreQualityByContext: StoreQualityByContextFromCmpH5.o
ccsh5tobam: CCSH5ToBam.o

CTESTS := \
ctest/alchemy.t  ctest/ccsh5tobam.t                  ctest/printTupleCountTable.t  ctest/sals.t      ctest/swmatcher.t \
ctest/bwt2sa.t   ctest/cmpH5StoreQualityByContext.t  ctest/sa2bwt.t                ctest/samodify.t


cramtests: ${EXE}
	cram -v --shell=/bin/bash ${CTESTS}

clean: 
	@rm -f ${EXE}
	@rm -f *.d *.o
