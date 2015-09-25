all:

SRCDIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))
-include ${CURDIR}/defines.mk
-include ${SRCDIR}/rules.mk
foo:
	echo $(realpath $(firstword $(MAKEFILE_LIST)))
	echo $(firstword $(MAKEFILE_LIST))
	echo $(MAKEFILE_LIST)
	echo ${SRCDIR}

CXXFLAGS += -O3 -g
CXXOPTS += \
		   -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive
GCXXFLAGS := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}
#INC_DIRS:=${LIBBLASR_INC} ${LIBPBIHDF_INC} ${LIBPBDATA_INC} ${PBBAM_INC} ${HTSLIB_INC} ${HDF5_INC} ${ZLIB_INC}
#LIB_DIRS:=${LIBBLASR_LIB} ${LIBPBIHDF_LIB} ${LIBPBDATA_LIB} ${PBBAM_LIB} ${HTSLIB_LIB} ${HDF5_LIB} ${ZLIB_LIB}
#LDLIBS := \
#	${LIBBLASR_LIBFLAGS} ${LIBPBIHDF_LIBFLAGS} ${LIBPBDATA_LIBFLAGS} \
#	${PBBAM_LIBFLAGS} ${HTSLIB_LIBFLAGS} ${HDF5_LIBFLAGS} ${ZLIB_LIBFLAGS} \
#	-ldl -lpthread

# HDF5 needs -ldl, but mobs does not pass it in.

SRCS := Blasr.cpp
OBJS := ${SRCS:.cpp=.o}
DEPS := ${SRCS:.cpp=.d}

LD_LIBRARY_PATH:=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
# Note: On macosx, this would be DYLD_LIBRARY_PATH.

vpath %.cpp ${SRCDIR}

init-submodule:
	${MAKE} update-submodule
	${MAKE} build-submodule

update-submodule:
	git submodule update --init

build-submodule:
	# DON'T use pbbam which is not on github.
	cd libcpp && NOPBBAM=true HDF5_LIB=${HDF5_LIB} HDF5_INC=${HDF5_INC} ./configure.py
	${MAKE} -C libcpp

submodule-clean:
	${RM} -r libcpp

# The rules above must be run separately.
all: blasr makeutils
#all: makeextrautils #This would require pbbam.
blasr: ${OBJS}
	${CXX} -o $@ ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${OBJS} ${LDFLAGS} ${LDLIBS}
	@echo LD_LIBRARY_PATH=${LD_LIBRARY_PATH}

makeutils:
	${MAKE} -C utils
makeextrautils:
	${MAKE} -C extrautils

CTESTS := \
ctest/affineAlign.t            ctest/bamOut.t    ctest/ccsH5.t            ctest/filtercriteria.t  ctest/m0-5.t             ctest/samNM.t \
ctest/aggressiveIntervalCut.t  ctest/bug25328.t  ctest/concordant.t       ctest/fofn.t            ctest/multipart.t        ctest/useccsallBestN1.t \
ctest/alignScore.t             ctest/bug25741.t  ctest/ecoli.t            ctest/hitpolicy.t       ctest/noSplitSubreads.t  ctest/useccsallLargeGenome.t\
ctest/bamIn.t                  ctest/bug25766.t  ctest/fastMaxInterval.t  ctest/holeNumbers.t     ctest/open_fail.t        ctest/verbose.t

SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

cramtests: blasr utils
	cram -v --shell=/bin/bash ${CTESTS}
	${MAKE} -C utils cramtests

cramfast: blasr utils
	cram -v --shell=/bin/bash $(filter-out ${SLOW_CTESTS},${CTESTS})
	${MAKE} -C utils cramfast

gtest: blasr
	# This requires the submodule to be configured with gtest.
	${MAKE} -C libcpp gtest

check: gtest cramtests

cleanall: cleanlib clean

# cleanlib is only for submodule users
cleanlib: libcpp/defines.mk
	${MAKE} -C libcpp clean

clean: 
	${RM} blasr ${OBJS} ${DEPS} blasr.d
	${MAKE} -C utils clean
	${MAKE} -C extrautils clean

-include ${DEPS}
