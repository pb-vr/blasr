#
# Definitions common to all make files.
#

HDF5INCLUDEDIR ?= ../../seymour/dist/common/include
HDF5LIBDIR     ?= ../../seymour/dist/common/lib/

INCLUDEDIRS = -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR)

HDF5LIB    ?= hdf5
HDF5LIBCPP ?= hdf5_cpp

GCCOPTS = -O3 -Wno-div-by-zero

HDF5LDFLAGS = -L$(HDF5LIBDIR) -l$(HDF5LIBCPP) -l$(HDF5LIB) -lz -lpthread -ldl

CC  ?= gcc
CXX ?= g++

CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS)
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  
CPP = $(CXX)

ifeq ($(shell $(CC) -dumpversion | awk -F '.' '$$1*100+$$2>404{print "yes"}'),yes)
    CPPOPTS += -fpermissive
endif

STATIC = 
# Set default value of STATIC to null. 
# In some operating systems, certain static libs (such as libz) may not be 
# support and blasr may fail to compile because of this.
# In order to compile blasr staticlly, please overwrite STATIC.
#ifneq ($(shell uname -s),Darwin)
#    STATIC   = -static
#endif
