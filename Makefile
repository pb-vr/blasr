SHELL=bash -e -E

.PHONY: all pblib makeutils gtest cramtests cramfast check clean cleanall cleanlib

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


all debug profile g: $(EXE) makeutils 

$(EXE): $(SRCS) $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

# DON'T use pbbam which is not on github.
PBINCROOT:=$(abspath libcpp)
pblib: $(PBINCROOT)/configure.py $(PBINCROOT)/makefile
	cd $(PBINCROOT) && NOPBBAM=true HDF5_LIB=${HDF5_LIB}/libhdf5.so ./configure.py
	cd $(PBINCROOT) && ${MAKE}

makeutils:
	make -C $(UTILS) $(MODE) 

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

cleanlib: $(PBINCROOT)/Makefile
	@make -C $(PBINCROOT) cleanall

clean: 
	@rm -f $(EXE) $(OBJS) $(DEPS) blasr.d
	@make -C $(UTILS) clean

-include $(DEPS)
