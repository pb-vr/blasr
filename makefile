all:

SRCDIR:=$(dir $(realpath $(firstword $(MAKEFILE_LIST))))
-include ${CURDIR}/defines.mk
-include ${SRCDIR}/rules.mk
foo:
	echo $(realpath $(firstword $(MAKEFILE_LIST)))
	echo $(firstword $(MAKEFILE_LIST))
	echo $(MAKEFILE_LIST)
	echo ${SRCDIR}

GET_SHA1 := $(shell git -C ${SRCDIR} describe --always --dirty='*')
CXXFLAGS += -O3 -g -DSHA1_7=\"${GET_SHA1}\"
CXXOPTS += \
		   -std=c++0x -pedantic \
           -Wall -Wextra -Wno-div-by-zero -Wno-overloaded-virtual \
           -MMD -MP
GCXXFLAGS := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
override CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}
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

override BLASR_PATH=${SRCDIR}/
export BLASR_PATH

override LD_LIBRARY_PATH:=${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}:${HDF5_LIB}:${HTSLIB_LIB}:${PBBAM_LIB}:${ZLIB_LIB}:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
# Note: On macosx, this would be DYLD_LIBRARY_PATH.

vpath %.cpp ${SRCDIR}

init-submodule:
	${MAKE} update-submodule
	${MAKE} configure-submodule
	${MAKE} build-submodule

update-submodule:
	git submodule update --init

configure-submodule:
	${MAKE} -f ${SRCDIR}/sub.mk configure-submodule

build-submodule:
	${MAKE} -C libcpp

distclean-submodule:
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

cramtests: blasr utils
	${MAKE} -f cram.mk cramtests
	${MAKE} -C utils cramtests

cramfast: blasr utils
	${MAKE} -f cram.mk cramfast
	${MAKE} -C utils cramfast

crammild: blasr utils
	${MAKE} -f cram.mk crammild
	${MAKE} -C utils crammild

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
