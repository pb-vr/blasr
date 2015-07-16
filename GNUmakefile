ROOT:=.
include q.mk
#include rules.mk

LDLIBS += -lpthread

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


all: blasr makeutils makeextrautils
blasr: ${OBJS}
	${CXX} -o $@ ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${OBJS} ${LDFLAGS} ${LDLIBS}

makeutils:
	${MAKE} -C utils
makeextrautils:
	${MAKE} -C extrautils

-include ${DEPS}
