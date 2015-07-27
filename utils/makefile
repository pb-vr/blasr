.PHONY=all cramtests

ROOT:=..
include ${ROOT}/defines.mk
include ${ROOT}/rules.mk

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}

EXE = loadPulses pls2fasta samtoh5 samtom4 samFilter toAfg sawriter sdpMatcher

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH


all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

loadPulses: LoadPulses.cpp
pls2fasta: PulseToFasta.cpp
samtoh5: SamToCmpH5.cpp
samtom4: SamToM4.cpp
samFilter: SamFilter.cpp
toAfg: ToAfg.cpp
sawriter: SAWriter.cpp
sdpMatcher: SDPMatcher.cpp

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/loadPulses.t ctest/pls2fasta.t

cramtests: ${EXE}
	cram -v --shell=/bin/bash ${CTESTS}

cramfast: ${EXE}
	cram -v --shell=/bin/bash ${filter-out ${SLOW_CTESTS}, ${CTESTS}}

clean: 
	@rm -f ${EXE}
	@rm -f *.d *.o
