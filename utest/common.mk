#
# Definitions common to all make files for unit test.
#

HDF5INCLUDEDIR ?= ../../seymour/dist/common/include
HDF5LIBDIR     ?= ../../seymour/dist/common/lib/
HDF5LIB    = hdf5
HDF5LIBCPP = hdf5_cpp

GTEST_DIR = /mnt/secondary-siv/third_party_source/gtest-1.6.0
# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h
GTEST_MAINA   = $(GTEST_DIR)/make/gtest_main.a  


INCLUDEDIRS = -I $(GTEST_DIR) -I $(GTEST_DIR)/include \
			  -I $(PBCPP_DIR)/common -I $(HDF5INCLUDEDIR)


GCCOPTS = -O3 -Wno-div-by-zero

CPPOPTS = $(GCCOPTS) $(INCLUDEDIRS)
CCOPTS  = $(GCCOPTS) $(INCLUDEDIRS)  

CPP = g++

LDFLAGS = -lpthread -L$(HDF5LIBDIR)



