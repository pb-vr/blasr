include defines.mk
include rules.mk

CXXOPTS := -O3 -g \
		   -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive
GCXXFLAGS := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
CXXFLAGS += ${CXXOPTS} ${GCXXFLAGS}
INC_DIRS:=${LIBBLASR_INCLUDE} ${LIBPBIHDF_INCLUDE} ${LIBPBDATA_INCLUDE} ${PBBAM_INCLUDE} ${HTSLIB_INCLUDE} ${HDF5_INCLUDE} ${ZLIB_INCLUDE}
LIB_DIRS:=${LIBBLASR_LIB} ${LIBPBIHDF_LIB} ${LIBPBDATA_LIB} ${PBBAM_LIB} ${HTSLIB_LIB} ${HDF5_LIB} ${ZLIB_LIB}
LDLIBS := \
	${LIBBLASR_LIBFLAGS} ${LIBPBIHDF_LIBFLAGS} ${LIBPBDATA_LIBFLAGS} \
	${PBBAM_LIBFLAGS} ${HTSLIB_LIBFLAGS} ${HDF5_LIBFLAGS} ${ZLIB_LIBFLAGS} \
	-ldl -lpthread
# HDF5 needs -ldl, but mobs does not pass it in.

SRCS := $(wildcard *.cpp)
OBJS := ${SRCS:.cpp=.o}
DEPS := ${SRCS:.cpp=.d}

all: blasr makeutils 

blasr: ${OBJS} ${PBLIB}
	${CXX} -o $@ ${CXXFLAGS} ${CPPFLAGS} -MF"${@:%=%.d}" ${OBJS} ${LDFLAGS} ${LDLIBS}

makeutils:
#	make -C utils ${MODE} 

-include ${DEPS}
