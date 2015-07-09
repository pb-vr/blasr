OS := $(shell uname)

# Determine where is PBINCROOT, either from github or PacBio SMRTAnalysis package.
PBINCROOT ?= $(realpath ../../lib/cpp/)
PREBUILT ?= $(realpath ../../../prebuilt.out)
THIRD_PARTY_PREFIX ?= $(realpath ..)

#####include $(PBINCROOT)/common.mk

LIBBLASR_INCLUDE  := $(PBINCROOT)/alignment
LIBPBIHDF_INCLUDE := $(PBINCROOT)/hdf
LIBPBDATA_INCLUDE := $(PBINCROOT)/pbdata
PBBAM_INCLUDE     := $(PBBAM)/include
HTSLIB_INCLUDE    := $(PBBAM)/../htslib/htslib
HDF5_INCLUDE      := $(HDF5_INC)
ZLIB_INCLUDE      := $(ZLIB_ROOT)/include

LIBBLASR_LIB  := $(PBINCROOT)/alignment
LIBPBIHDF_LIB := $(PBINCROOT)/hdf
LIBPBDATA_LIB := $(PBINCROOT)/pbdata
PBBAM_LIB     := $(PBBAM)/lib
HTSLIB_LIB    := $(PBBAM)/../htslib
HDF5_LIB      := $(HDF5_LIB)
ZLIB_LIB      := $(ZLIB_ROOT)/lib


LIBPBIHDF_CPP_LIBFLAG := -lhdf5_cpp
LIBPBIHDF_LIBFLAG     := -lhdf5

LIBBLASR_LIBFLAGS  := -lblasr
LIBPBIHDF_LIBFLAGS := -lpbihdf
LIBPBDATA_LIBFLAGS := -lpbdata
PBBAM_LIBFLAGS     := -lpbbam
HTSLIB_LIBFLAGS    := -lhts
HDF5_LIBFLAGS      := $(LIBHDF5_CPP_LIBFLAG) $(LIBHDF5_LIBFLAG)
ZLIB_LIBFLAGS      := -lz


INCDIRS = -I$(LIBBLASR_INCLUDE) \
          -I$(LIBPBIHDF_INCLUDE) \
          -I$(LIBPBDATA_INCLUDE) \
          -I$(HDF5_INCLUDE)

LIBDIRS = -L$(LIBBLASR_LIB) \
          -L$(LIBPBIHDF_LIB) \
          -L$(LIBPBDATA_LIB) \
          -L$(HDF5_LIB) 

ifneq ($(ZLIB_ROOT), notfound)
        # NOTE: The zlib include directory is not needed since we do not
        #       refer to any zlib include files directly in the blasr source
        #       code.   But we do need to link to the zlib library (static
        #       or shared), since htslib depends on it.
        #          INCDIRS += -I$(ZLIB_INCLUDE)
	LIBDIRS += -L$(ZLIB_LIB)
endif

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive

DEFAULTCXXFLAG := -O3
DEBUGCXXFLAG := -g -ggdb -fno-inline
PROFILECXXFLAG := -Os -pg 
GCXXFLAG := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 

LIBS := $(LIBBLASR_LIBFLAGS) $(LIBPBIHDF_LIBFLAGS) $(LIBPBDATA_LIBFLAGS)
PBBAMLIBS :=  $(PBBAM_LIBFLAGS)
ifneq ($(SHARED_LIB),)
    PBBAMLIBS += $(HTSLIB_LIBFLAGS)
else
    ifneq ($(wildcard "$(HTSLIB_LIB)/libhts.a"),"")
        PBBAMLIBS += $(HTSLIB_LIB)/libhts.a
    else
         PBBAMLIBS += $(HTSLIB_LIBFLAGS)
    endif
endif

ifneq ($(SHARED_LIB),)
    LIBS += $(HDF5_LIBFLAGS)
else
    ifneq ($(wildcard "$(HDF5_LIB)/libhdf5_cpp.a"),"")
        LIBS += $(HDF5_LIB)/libhdf5_cpp.a $(HDF5_LIB)/libhdf5.a
    else
        LIBS += $(HDF5_LIBFLAGS)
    endif
endif

ifeq ($(origin nopbbam), undefined)
	INCDIRS += -I$(PBBAM_INCLUDE) \
			   -I$(HTSLIB_INCLUDE) \
			   -I$(BOOST_INCLUDE)


	LIBDIRS += -L$(PBBAM_LIB) \
			   -L$(HTSLIB_LIB)

	LIBS += $(PBBAMLIBS) 
endif

LIBS += $(ZLIB_LIBFLAGS)
LIBS += -lpthread

# NOTE: libdl is not needed when we are linking with shared libraries, but
#       is is needed for static linking, since the hdf5 library depends on it
ifeq ($(SHARED_LIB),)
    LIBS += -ldl
endif

ifneq ($(SHARED_LIB),)
    STATIC :=
else
    ifneq ($(OS), Darwin)
	STATIC := -static
    else
	STATIC :=
    endif
endif

ifneq ($(OS), Darwin)
	LIBS += -lrt
else
	LIBS += -lsz
endif

# -lhdf5, -lhdf5_cpp, -lz required for HDF5
# -lpbbam -lhts for BAM
# -lpthread for multi-threading
# -lrt for clock_gettime
# -ldl for dlopen dlclose 

GLIBS = -Wl --eh-frame-hdr -fno-builtin-malloc -L$(HOME)/lib -ltcmalloc -lunwind -lprofiler $(LIBS)

PBLIB :=
ifeq ($(NO_SUBMAKES),)
    PBLIB := pblib
endif

