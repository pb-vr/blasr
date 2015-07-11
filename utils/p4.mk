# This p4.mk is only used to make 
#   //depot/software/smrtanalysis/bioinformatics/tools/blasr
# on PacBio internal p4 build.

SHELL=/bin/bash -e -E

.PHONY=all cramtests

PBINCROOT ?= $(realpath ../../../lib/cpp/)
PREBUILT ?= $(realpath ../../../../prebuilt.out)
THIRD_PARTY_PREFIX := $(realpath ../..)

include ../p4.common.mk

all : CXXFLAGS ?= $(DEFAULTCXXFLAG)
debug : CXXFLAGS ?= $(DEBUGCXXFLAG)
profile : CXXFLAGS ?= $(PROFILECXXFLAG)
g: CXXFLAGS += $(GCXXFLAG)
g: LIBS = $(GLIBS)

EXE = loadPulses pls2fasta samtoh5 samtom4 samFilter toAfg sawriter sdpMatcher

all debug profile g: $(EXE)

loadPulses: LoadPulses.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

pls2fasta: PulseToFasta.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

samtoh5: SamToCmpH5.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

samtom4: SamToM4.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

samFilter: SamFilter.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

toAfg: ToAfg.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

sawriter: SAWriter.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

sdpMatcher: SDPMatcher.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

pblib: 
	make -C ../ -f p4.mk mkliba

CTESTS := $(wildcard ctest/*.t)
SLOW_CTESTS := ctest/loadPulses.t ctest/pls2fasta.t

cramtests: $(EXE)
	cram -v --shell=/bin/bash $(CTESTS)

cramfast: $(EXE)
	cram -v --shell=/bin/bash $(filter-out $(SLOW_CTESTS), $(CTESTS))

clean: 
	@rm -f $(EXE)
	@rm -f *.d *.o

-include $(DEPS)
