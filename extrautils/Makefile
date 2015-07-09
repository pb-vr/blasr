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

exe = sa2bwt bwt2sa alchemy excrep evolve bsdb simpleShredder swMatcher \
  samodify sals printTupleCountTable cmpH5StoreQualityByContext ccsh5tobam

all debug profile g: $(exe)

sa2bwt: SuffixArrayToBWT.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

bwt2sa: BwtToSuffixArray.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

alchemy: BasH5Simulator.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

excrep: ExciseRepeats.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

evolve: Evolve.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

bsdb: BuildSequenceDB.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

simpleShredder: SimpleShredder.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

swMatcher: SWMatcher.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

samodify: SAModify.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

sals: SALS.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

printTupleCountTable: PrintTupleCountTable.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

cmpH5StoreQualityByContext: StoreQualityByContextFromCmpH5.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

ccsh5tobam: CCSH5ToBam.cpp $(PBLIB)
	$(CXX) $(CXXOPTS) $(CXXFLAGS) $(INCDIRS) -MF"$(@:%=%.d)" $(STATIC) -o $@ $< $(LIBDIRS) $(LIBS)

pblib: $(PBINCROOT)/Makefile
	make -C $(PBINCROOT)

.INTERMEDIATE: $(OBJS)

CTESTS = $(wildcard ctest/*.t)
cramtests: $(EXE)
	cram --shell=/bin/bash $(CTESTS)

clean: 
	@rm -f $(exe)
	@rm -f $(OBJS) $(DEPS)
	@rm *.o *.d

-include $(DEPS)
