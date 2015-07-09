SHELL=/bin/bash -e -E

.PHONY: all pblib makeutils gtest cramtests cramfast check clean cleanall cleanlib

nopbbam=true

GIT_BLASR_LIBPATH = libcpp
PB_BLASR_LIBPATH = ../../lib/cpp
# Determine where is PBINCROOT, either from github or PacBio SMRTAnalysis package.
PBINCROOT ?= $(shell cd $(GIT_BLASR_LIBPATH) 2>/dev/null && pwd || echo -n notfound)
ifeq ($(PBINCROOT), notfound)
	PBINCROOT := $(shell cd $(PB_BLASR_LIBPATH) 2>/dev/null && pwd || echo -n notfound)
	ifeq ($(PBINCROOT), notfound)
		$(error please check your blasr lib exists.)
	endif
endif

# common.mk contains the configuration for this build setup
GIT_COMMON_MK = blasr_git_common.mk
ifneq ($(shell ls $(GIT_COMMON_MK) 2>/dev/null || echo -n notfound), notfound)
include $(GIT_COMMON_MK)
endif

include common.mk

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)
EXE = blasr
UTILS = utils

all : CXXFLAGS ?= $(DEFAULTCXXFLAG)
debug : CXXFLAGS ?= $(DEBUGCXXFLAG)
profile : CXXFLAGS ?= $(PROFILECXXFLAG)
g: CXXFLAGS += $(GCXXFLAG)
g: LIBS = $(GLIBS)

all: MODE = 
debug: MODE=debug
profile: MODE = profile 
g: MODE = g

all debug profile g: $(EXE) makeutils

$(EXE): $(SRCS) $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

# DON'T use pbbam which is not on github.
pblib: $(PBINCROOT)/configure.py $(PBINCROOT)/makefile
	cd $(PBINCROOT) && NOPBBAM=true HDF5_LIB=${HDF5_LIB}/libhdf5.so ./configure.py
	cd $(PBINCROOT) && make -j

makeutils:
	export PBINCROOT=$(PBINCROOT) && export nopbbam=true && export COMMON_NO_THIRD_PARTY_REQD=true && export HDF5_LIB=$(HDF5_LIB) && export HDF5_INC=$(HDF5_INC) && make -C $(UTILS) $(MODE) 

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

cramtests: $(EXE) $(UTILS)
	cram -v --shell=/bin/bash $(CTESTS)
	make -C $(UTILS) cramtests

cramfast: $(EXE) $(UTILS)
	cram -v --shell=/bin/bash $(filter-out $(SLOW_CTESTS),$(CTESTS))
	make -C $(UTILS) cramfast

gtest: $(EXE)
	make -C $(PBINCROOT) gtest

check: gtest cramtests

cleanall: cleanlib clean

cleanlib: $(PBINCROOT)/makefile
	@COMMON_NO_THIRD_PARTY_REQD=true make -C $(PBINCROOT) cleanall

clean: 
	@rm -f $(EXE) $(OBJS) $(DEPS) blasr.d
	@COMMON_NO_THIRD_PARTY_REQD=true PBINCROOT=$(PBINCROOT) make -C $(UTILS) clean

-include $(DEPS) 
