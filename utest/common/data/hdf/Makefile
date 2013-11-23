.PHONY : all gtest clean
ROOT_DIR = ../../..
PBCPP_DIR = $(ROOT_DIR)/..
include $(ROOT_DIR)/common.mk

HDF5LIBDIR ?= $(ROOT_DIR)/$(H5DIR)/lib/
HDF5INCLUDEDIR ?= $(ROOT_DIR)/$(H5DIR)/include/

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECS = $(SOURCES:.cpp=)

CCOPTS += -I$(HDF5INCLUDEDIR)
LDFLAGS += -L$(HDF5LIBDIR) -l$(HDF5LIB) -l$(HDF5LIBCPP) -lz 


all: $(SOURCES) $(OBJECTS) $(EXECS)

$(EXECS): % : %.o
	$(CXX) $< $(GTEST_MAINA) -o $@ $(LDFLAGS)

.cpp.o: 
	$(CXX) $(CCOPTS) -c $< -o $@


gtest: $(EXECS)
	@(for exe in $(EXECS); do echo "Running $$exe"; ./$$exe --gtest_output=xml:$(ROOT_DIR)/$$exe.xml; done)

clean: 
	@rm -f $(OBJECTS) $(EXECS)

