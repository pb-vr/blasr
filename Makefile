SHELL=bash

.PHONY: all pblib gtest cramtests check clean cleanall

include common.mk

SRCS := $(wildcard *.cpp)
DEPS := $(SRCS:.cpp=.d)
EXE = blasr

all : CXXFLAGS ?= $(DEFAULTCXXFLAG)
debug : CXXFLAGS ?= $(DEBUGCXXFLAG)
profile : CXXFLAGS ?= $(PROFILECXXFLAG)
g: CXXFLAGS += $(GCXXFLAG)
g: LIBS = $(GLIBS)

UTILS = utils

all debug profile g: $(EXE) $(UTILS)
	make -C $(UTILS)

$(EXE): $(SRCS) $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

pblib: $(PBINCROOT)/Makefile
	make -C $(PBINCROOT)

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

cramtests: $(EXE)
	cram -v --shell=/bin/bash $(CTESTS)
	make -C $(UTILS) cramtests

cramfast: $(EXE)
	cram -v --shell=/bin/bash $(filter-out $(SLOW_CTESTS),$(CTESTS))
	make -C $(UTILS) cramfast

gtest: $(EXE)
	make -C $(PBINCROOT) gtest

check: gtest cramtests

clean: 
	@rm -f $(EXE)
	@rm -f $(DEPS)
	@rm -f blasr.d
	@make -C $(UTILS) clean

cleanall: clean $(PBINCROOT)/Makefile
	make -C $(PBINCROOT) cleanall

-include $(DEPS)
