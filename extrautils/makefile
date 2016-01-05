.PHONY=all cramtests

SRCDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
-include ${CURDIR}/../defines.mk
include ${SRCDIR}/../rules.mk

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wextra \
           -MMD -MP
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}

EXE = sa2bwt bwt2sa excrep evolve bsdb simpleShredder swMatcher \
  samodify sals printTupleCountTable cmpH5StoreQualityByContext

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH

vpath %.cpp ${SRCDIR}

all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

sa2bwt: SuffixArrayToBWT.o
bwt2sa: BwtToSuffixArray.o
excrep: ExciseRepeats.o
evolve: Evolve.o
bsdb: BuildSequenceDB.o
simpleShredder: SimpleShredder.o
swMatcher: SWMatcher.o
samodify: SAModify.o
sals: SALS.o
printTupleCountTable: PrintTupleCountTable.o
cmpH5StoreQualityByContext: StoreQualityByContextFromCmpH5.o

CTESTS := \
ctest/printTupleCountTable.t  ctest/sals.t      ctest/swmatcher.t \
ctest/bwt2sa.t   ctest/cmpH5StoreQualityByContext.t  ctest/sa2bwt.t                ctest/samodify.t


cramtests: ${EXE}
	cram -v --shell=/bin/bash ${CTESTS}

clean: 
	@rm -f ${EXE}
	@rm -f *.d *.o
