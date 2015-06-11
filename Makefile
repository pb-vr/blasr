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
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

ifeq ($(origin nopbbam), undefined)
pblib: $(PBINCROOT)/Makefile
	@echo building pblib with pbbam
	make -C $(PBINCROOT)
else
pblib: $(PBINCROOT)/Makefile
	@echo building pblib without pbbam
	nopbbam=true make -C $(PBINCROOT)
endif

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
	@rm -f $(EXE) $(OBJS) $(DEPS)
	@make -C $(UTILS) clean

-include $(DEPS)
