PBINCROOT ?= ../../lib/cpp
HDFINC ?= ../../../assembly/seymour/dist/common/include
HDFLIB ?= ../../../assembly/seymour/dist/common/lib

INCDIRS = -I$(PBINCROOT)/alignment \
		  -I$(PBINCROOT)/pbdata \
		  -I$(PBINCROOT)/hdf \
		  -I$(HDFINC) 

LIBDIRS = -L$(PBINCROOT)/alignment \
		  -L$(PBINCROOT)/pbdata \
		  -L$(PBINCROOT)/hdf \
		  -L$(HDFINC) \
		  -L$(HDFLIB)

CXXFLAGS := -std=c++0x -Wall -Wuninitialized -Wno-div-by-zero \
			-pedantic -c -fmessage-length=0 -MMD -MP -w -fpermissive

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(SRCS:.cpp=.d)
LIBS := -lblasr -lpbdata -lpbihdf -lhdf5_cpp -lhdf5 -lz -lpthread -lrt -ldl
# -lhdf5, -lhdf5_cpp, -lz required for HDF5
# -lpthread for multi-threading
# -lrt for clock_gettime
# -ldl for dlopen dlclose 


all : OPTIMIZE = -O3

debug : OPTIMIZE = -g -ggdb -fno-inline

profile : OPTIMIZE = -Os -pg 

g: OPTIMIZE = -O3
g: G_OPTIMIZE = -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer 
g: G_LIBS = -Wl --eh-frame-hdr -fno-builtin-malloc -L$(HOME)/lib -ltcmalloc -lunwind -lprofiler $(LIBS)


exe=blasr
all debug profile g: $(exe) 

$(exe): Blasr.o
	$(CXX) $(LIBDIRS) $(OPTIMIZE) $(G_OPTIMIZE) -static -o $@ $^ $(LIBS) $(G_LIBS)

Blasr.o: Blasr.cpp
	$(CXX) $(CXXFLAGS) $(OPTIMIZE) $(G_OPTIMIZE) $(INCDIRS) -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.o) $(@:%.o=%.d)" -o $@ $<

.INTERMEDIATE: $(OBJS)

cramtests:
	cram --shell=/bin/bash ctest/*.t

clean: 
	rm -f $(exe)
	rm -f $(OBJS) $(DEPS)

-include $(DEPS)
