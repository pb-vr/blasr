SHELL=/bin/bash -e -E

.PHONY=all cramtests

PBINCROOT ?= $(realpath ../../../lib/cpp/)
PREBUILT ?= $(realpath ../../../../prebuilt.out)
THIRD_PARTY_PREFIX := $(realpath ../..)

include ../common.mk

all : CXXFLAGS ?= $(DEFAULTCXXFLAG)

debug : CXXFLAGS ?= $(DEBUGCXXFLAG)

profile : CXXFLAGS ?= $(PROFILECXXFLAG)

g: CXXFLAGS += $(GCXXFLAG)
g: LIBS = $(GLIBS)

EXE = loadPulses pls2fasta samtoh5 samtom4 samFilter toAfg sawriter sdpmatcher

all debug profile g: $(EXE)

loadPulses: LoadPulses.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ LoadPulses.cpp $(LIBDIRS) $(LIBS)

pls2fasta: PulseToFasta.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ PulseToFasta.cpp $(LIBDIRS) $(LIBS)

samtoh5: SamToCmpH5.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ SamToCmpH5.cpp $(LIBDIRS) $(LIBS)

samtom4: SamToM4.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ SamToM4.cpp $(LIBDIRS) $(LIBS)

samFilter: SamFilter.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ SamFilter.cpp $(LIBDIRS) $(LIBS)

toAfg: ToAfg.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ ToAfg.cpp $(LIBDIRS) $(LIBS)

sawriter: SAWriter.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ SAWriter.cpp $(LIBDIRS) $(LIBS)

sdpmatcher: SDPMatcher.cpp $(PBLIB)
	$(CXX_pp) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ SDPMatcher.cpp $(LIBDIRS) $(LIBS)

pblib: $(PBINCROOT)/Makefile
	make -C $(PBINCROOT)

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
