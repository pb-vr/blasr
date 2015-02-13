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
          -I$(HDF5_ROOT)/src \
          -I$(HDF5_ROOT)/c++/src

LIBDIRS = -L$(PBINCROOT)/alignment \
          -L$(PBINCROOT)/pbdata \
          -L$(PBINCROOT)/hdf \
          -L$(HDF5_ROOT)/src/.libs \
          -L$(HDF5_ROOT)/c++/src/.libs

ifneq ($(ZLIB_ROOT), notfound)
	INCDIRS += -I$(ZLIB_ROOT)/include
	LIBDIRS += -L$(ZLIB_ROOT)/lib
endif

CXXOPTS := -std=c++0x -pedantic \
           -Wall -Wuninitialized -Wno-div-by-zero \
           -MMD -MP -w -fpermissive

SRCS := $(wildcard *.cpp)
DEPS := $(SRCS:.cpp=.d)
LIBS := -lblasr -lpbdata -lpbihdf -lhdf5_cpp -lhdf5 -lz -lpthread -ldl
EXE  := blasr

ifneq ($(OS), Darwin)
	LIBS += -lrt
	STATIC := -static
else
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

cramtests: $(EXE)
	cram -v --shell=/bin/bash ctest/*.t

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
