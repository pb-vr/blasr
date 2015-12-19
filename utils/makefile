.PHONY=all cramtests

SRCDIR:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
-include ${CURDIR}/../defines.mk
include ${SRCDIR}/../rules.mk

CXXFLAGS ?= -O3 -g
CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wextra \
           -MMD -MP -w
override CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}

EXE = loadPulses pls2fasta samtoh5 samtom4 samFilter toAfg sawriter sdpMatcher

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH

vpath %.cpp ${SRCDIR}

all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

loadPulses: LoadPulses.o
pls2fasta: PulseToFasta.o
samtoh5: SamToCmpH5.o
samtom4: SamToM4.o
samFilter: SamFilter.o
toAfg: ToAfg.o
sawriter: SAWriter.o
sdpMatcher: SDPMatcher.o

CTESTS := \
  ctest/loadPulses_fast.t  ctest/pls2fasta_fast.t  ctest/samFilter.t  ctest/samtom4.t   ctest/sdpMatcher.t \
  ctest/loadPulses.t       ctest/pls2fasta.t       ctest/samtoh5.t    ctest/sawriter.t  ctest/toAfg.t

SLOW_CTESTS := ctest/loadPulses.t ctest/pls2fasta.t

cramtests: ${EXE}
	cram -v --shell=/bin/bash ${CTESTS}

cramfast: ${EXE}
	cram -v --shell=/bin/bash ${filter-out ${SLOW_CTESTS}, ${CTESTS}}

clean: 
	@rm -f ${EXE}
	@rm -f *.d *.o
