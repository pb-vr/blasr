all:

ROOT:=.
include defines.mk
include rules.mk

CXXFLAGS += -O3 -g
CXXOPTS += \
		   -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive
GCXXFLAGS := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}
#INC_DIRS:=${LIBBLASR_INCLUDE} ${LIBPBIHDF_INCLUDE} ${LIBPBDATA_INCLUDE} ${PBBAM_INCLUDE} ${HTSLIB_INCLUDE} ${HDF5_INCLUDE} ${ZLIB_INCLUDE}
#LIB_DIRS:=${LIBBLASR_LIB} ${LIBPBIHDF_LIB} ${LIBPBDATA_LIB} ${PBBAM_LIB} ${HTSLIB_LIB} ${HDF5_LIB} ${ZLIB_LIB}
#LDLIBS := \
#	${LIBBLASR_LIBFLAGS} ${LIBPBIHDF_LIBFLAGS} ${LIBPBDATA_LIBFLAGS} \
#	${PBBAM_LIBFLAGS} ${HTSLIB_LIBFLAGS} ${HDF5_LIBFLAGS} ${ZLIB_LIBFLAGS} \
#	-ldl -lpthread

# HDF5 needs -ldl, but mobs does not pass it in.

SRCS := $(wildcard *.cpp)
OBJS := ${SRCS:.cpp=.o}
DEPS := ${SRCS:.cpp=.d}

LD_LIBRARY_PATH=${HDF5_LIB}:${LIBBLASR_LIB}:${LIBPBIHDF_LIB}:${LIBPBDATA_LIB}
export LD_LIBRARY_PATH


init-submodule:
	${MAKE} update-submodule
	${MAKE} build-submodule

update-submodule:
	git submodule update --init

build-submodule:
	# DON'T use pbbam which is not on github.
	cd libcpp && NOPBBAM=true HDF5_LIB=${HDF5_LIB}/libhdf5.so HDF5_INCLUDE=${HDF5_INCLUDE} ./configure.py
	${MAKE} -C libcpp

submodule-clean:
	${RM} -r blasr_libcpp

# The rules above must be run separately.
all: blasr makeutils makeextrautils
blasr: ${OBJS}
	${CXX} -o $@ ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${OBJS} ${LDFLAGS} ${LDLIBS}
	echo ${LD_LIBRARY_PATH}

makeutils:
	${MAKE} -C utils
makeextrautils:
	${MAKE} -C extrautils

-include ${DEPS}
