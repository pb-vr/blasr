OS := $(shell uname)

# Determine where is PBINCROOT, either from github or PacBio SMRTAnalysis package.
PBINCROOT ?= $(realpath ../../lib/cpp/)
PREBUILT ?= $(realpath ../../../prebuilt.out)
THIRD_PARTY_PREFIX ?= $(realpath ..)

include $(PBINCROOT)/common.mk

LIBBLASR_INCLUDE  := $(PBINCROOT)/alignment
LIBPBIHDF_INCLUDE := $(PBINCROOT)/hdf
LIBPBDATA_INCLUDE := $(PBINCROOT)/pbdata
PBBAM_INCLUDE     := $(PBBAM)/include
HTSLIB_INCLUDE    := $(PBBAM)/third-party/htslib

LIBBLASR_LIB  := $(PBINCROOT)/alignment
LIBPBIHDF_LIB := $(PBINCROOT)/hdf
LIBPBDATA_LIB := $(PBINCROOT)/pbdata
PBBAM_LIB     := $(PBBAM)/lib
HTSLIB_LIB    := $(PBBAM)/third-party/htslib

INCDIRS = -I$(LIBBLASR_INCLUDE) \
          -I$(LIBPBIHDF_INCLUDE) \
          -I$(LIBPBDATA_INCLUDE) \
          -I$(HDF5_INC) \
          -I$(PBBAM_INCLUDE) \
          -I$(HTSLIB_INCLUDE) \
          -I$(BOOST_INCLUDE)

LIBDIRS = -L$(LIBBLASR_LIB) \
          -L$(LIBPBIHDF_LIB) \
          -L$(LIBPBDATA_LIB) \
          -L$(HDF5_LIB) \
          -L$(PBBAM_LIB) \
          -L$(HTSLIB_LIB)

ifneq ($(ZLIB_ROOT), notfound)
	INCDIRS += -I$(ZLIB_ROOT)/include
	LIBDIRS += -L$(ZLIB_ROOT)/lib
endif

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive

DEFAULTCXXFLAG := -O3
DEBUGCXXFLAG := -g -ggdb -fno-inline
PROFILECXXFLAG := -Os -pg 
GCXXFLAG := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 

LIBS := -lblasr -lpbihdf -lpbdata -lpbbam 
ifneq ($(wildcard "$(HTSLIB_LIB)/libhts.a"),"")
	LIBS += $(HTSLIB_LIB)/libhts.a
else
	LIBS += -lhts
endif
ifneq ($(wildcard "$(HDF5_LIB)/libhdf5_cpp.a"),"")
    LIBS += $(HDF5_LIB)/libhdf5_cpp.a $(HDF5_LIB)/libhdf5.a -lz -lpthread -ldl
else
    LIBS += -lhdf5_cpp -lhdf5 -lz -lpthread -ldl
endif

ifneq ($(OS), Darwin)
	LIBS += -lrt
	STATIC := -static
else
	LIBS += -lsz
	STATIC :=
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

