SHELL=bash

.PHONY: all pblib gtest cramtests check clean cleanall

OS := $(shell uname)

GIT_BLASR_LIBPATH = ../lib
PB_BLASR_LIBPATH = ../../lib/cpp

# Determine where is PBINCROOT, either from github or PacBio SMRTAnalysis package.
PBINCROOT ?= $(shell cd $(GIT_BLASR_LIBPATH) 2>/dev/null && pwd || echo -n notfound)
ifeq ($(PBINCROOT), notfound)
	PBINCROOT := $(shell cd $(PB_BLASR_LIBPATH) 2>/dev/null && pwd || echo -n notfound)
	ifeq ($(PBINCROOT), notfound)
		$(error please check your blasr lib exists.)
	endif
endif

PREBUILT ?= ../../../prebuilt.out
THIRD_PARTY_PREFIX := ..

include $(PBINCROOT)/common.mk

INCDIRS = -I$(PBINCROOT)/alignment \
          -I$(PBINCROOT)/pbdata \
          -I$(PBINCROOT)/hdf \
          -I$(HDF5_INC) \
          -I$(PBBAM)/include \
          -I$(PBBAM)/third-party/htslib \
          -I$(PREBUILT)/boost/boost_1_55_0

LIBDIRS = -L$(PBINCROOT)/alignment \
          -L$(PBINCROOT)/pbdata \
          -L$(PBINCROOT)/hdf \
          -L$(HDF5_LIB) \
          -L$(PBBAM)/lib \
          -L$(PBBAM)/third-party/htslib

ifneq ($(ZLIB_ROOT), notfound)
	INCDIRS += -I$(ZLIB_ROOT)/include
	LIBDIRS += -L$(ZLIB_ROOT)/lib
endif

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive

SRCS := $(wildcard *.cpp)
DEPS := $(SRCS:.cpp=.d)
ifneq ($(wildcard "$(HDF5_LIB)/libhdf5_cpp.a"),"")
    LIBS := -lblasr -lpbdata -lpbihdf -lpbbam -lhts $(HDF5_LIB)/libhdf5_cpp.a $(HDF5_LIB)/libhdf5.a -lz -lpthread -ldl
else
    LIBS := -lblasr -lpbdata -lpbihdf -lpbbam -lhts -lhdf5_cpp -lhdf5 -lz -lpthread -ldl
endif
EXE  := blasr

ifneq ($(OS), Darwin)
	LIBS += -lrt
	STATIC := -static
else
	LIBS += -lsz
	STATIC :=
endif

# -lhdf5, -lhdf5_cpp, -lz required for HDF5
# -lpthread for multi-threading
# -lrt for clock_gettime
# -ldl for dlopen dlclose 


all : CXXFLAGS ?= -O3

debug : CXXFLAGS ?= -g -ggdb -fno-inline

profile : CXXFLAGS ?= -Os -pg 

g: CXXFLAGS += -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
g: LIBS += -Wl --eh-frame-hdr -fno-builtin-malloc -L$(HOME)/lib -ltcmalloc -lunwind -lprofiler $(LIBS)

all debug profile g: $(EXE) 

$(EXE): $(SRCS) pblib
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

pblib: $(PBINCROOT)/Makefile
	make -C $(PBINCROOT)

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

cramtests: $(EXE)
	cram -v --shell=/bin/bash $(CTESTS)

cramfast: $(EXE)
	cram -v --shell=/bin/bash $(filter-out $(SLOW_CTESTS),$(CTESTS))

gtest: $(EXE)
	make -C $(PBINCROOT) gtest

check: gtest cramtests

clean: 
	@rm -f $(EXE)
	@rm -f $(DEPS)
	@rm -f blasr.d

cleanall: clean $(PBINCROOT)/Makefile
	make -C $(PBINCROOT) cleanall

-include $(DEPS)
