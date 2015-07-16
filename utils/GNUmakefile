SHELL=/bin/bash -e -E

.PHONY=all cramtests

ROOT:=..
include q.mk

INCDIRS := \
	${LIBBLASR_INCLUDE} \
	${LIBPBDATA_INCLUDE} \
	${LIBPBIHDF_INCLUDE} \
	${PBBAM_INCLUDE} \
	${HDF5_INCLUDE} \
	${HTSLIB_INCLUDE} \
	${BOOST_INCLUDE}
LIBDIRS := \
	${LIBBLASR_LIB} \
	${LIBPBDATA_LIB} \
	${LIBPBIHDF_LIB} \
	${PBBAM_LIB} \
	${HDF5_LIB} \
	${HTSLIB_LIB} \
	${ZLIB_LIB}
LDLIBS+= \
	${LIBBLASR_LIBFLAGS} \
	${LIBPBDATA_LIBFLAGS} \
	${LIBPBIHDF_LIBFLAGS} \
	${PBBAM_LIBFLAGS} \
	${HDF5_LIBFLAGS} \
	${HTSLIB_LIBFLAGS} \
	${ZLIB_LIBFLAGS} \
	${PTHREAD_LIBFLAGS} \
	${DL_LIBFLAGS}
LDLIBS+=-lrt
CPPFLAGS+=$(patsubst %,-I%,${INCDIRS})
LDFLAGS+=$(patsubst %,-L %,${LIBDIRS})

ALL_CXXFLAGS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive \
		   ${CXXFLAGS}

EXE = loadPulses pls2fasta samtoh5 samtom4 samFilter toAfg sawriter sdpMatcher

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH


all: ${EXE}

${EXE}:
	${CXX} -o $@ $< ${ALL_CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${STATIC} ${LDFLAGS} ${LDLIBS}

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

-include ${DEPS}
