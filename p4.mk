SHELL=bash -e -E -vx

.PHONY: all pblib makeutils gtest cramtests cramfast check clean wipe

include p4.common.mk

BLASR_LIBCPP := blasr_libcpp
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
MAKE := make

all debug profile g: $(EXE) makeutils 

$(EXE): $(SRCS) mkliba
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $(SRCS) $(LIBDIRS) $(LIBS)

makeutils: mkliba
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(UTILS) $(MODE) -f p4.mk

fetch: $(PBINCROOT)
	@echo Fetching $(BLASR_LIBCPP)
	@echo "SHELL=bash -evx" > $(FETCHMK) && grep "^REF:=" $(SUBMK) >> $(FETCHMK) && grep "^NAME:=" $(SUBMK) >> $(FETCHMK) && echo "REPO_DIR:=$(BLASR_LIBCPP)" >> $(FETCHMK) && cat  $(FETCHTMP) >> $(FETCHMK)
	make -f $(FETCHMK)

mkliba: fetch mkpbbama
	$(BLASR_LIBCPP)/configure.py BOOST_INCLUDE=$(BOOST_INCLUDE) HTSLIB_INCLUDE=$(HTSLIB_INCLUDE) PBBAM_INCLUDE=$(PBBAM_INCLUDE) HTSLIB_LIB=$(HTSLIB_LIB) PBBAM_LIB=$(PBBAM_LIB)
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(BLASR_LIBCPP)/pbdata/ libpbdata.a
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(BLASR_LIBCPP)/hdf/ libpbihdf.a
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(BLASR_LIBCPP)/alignment/ libblasr.a

mklibso: fetch mkpbbamso
	$(BLASR_LIBCPP)/configure.py SHARED_LIB=true BOOST_INCLUDE=$(BOOST_INCLUDE) HTSLIB_INCLUDE=$(HTSLIB_INCLUDE) PBBAM_INCLUDE=$(PBBAM_INCLUDE) HTSLIB_LIB=$(HTSLIB_SO) PBBAM_LIB=$(PBBAM_SO)
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(BLASR_LIBCPP) all

mkpbbama: 
	ls $(PBBAM_A) 2>/dev/null 1>/dev/null || (cd $(PBBAM) && rm -f builda && mkdir builda && cd builda && cmake ../ && make)

mkpbbamso:
	ls $(PBBAM_SO) 2>/dev/null 1>/dev/null || cd $(PBBAM) && rm -f buildso && mkdir -p buildso && cd buildso && cmake -DPacBioBAM_build_shared=ON ../ && make 

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/bug25328.t ctest/useccsallLargeGenome.t

cramtests: $(EXE) $(UTILS) 
	cram -v --shell=/bin/bash $(CTESTS)
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(UTILS) cramtests

cramfast: $(EXE) $(UTILS)
	cram -v --shell=/bin/bash $(filter-out $(SLOW_CTESTS),$(CTESTS))
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(UTILS) cramfast

gtest: $(EXE)
	export CXXFLAGS="$(CXXFLAGS)" && make -C $(PBINCROOT) gtest

check: gtest cramtests

wipe: clean
	@echo WIPING OUT $(BLASR_LIBCPP)...
	@rm -rf $(BLASR_LIBCPP)

clean: 
	@rm -f $(EXE) $(OBJS) $(DEPS) blasr.d 
	@make -C $(UTILS) -f p4.mk clean
	@make -C $(BLASR_LIBCPP) cleanall
	@rm -f $(FETCHMK)

-include $(DEPS)
