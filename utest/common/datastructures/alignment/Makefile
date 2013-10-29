ROOT_DIR = ../../..
PBCPP_DIR = $(ROOT_DIR)/..
include $(ROOT_DIR)/common.mk

LDFLAGS += -l$(HDF5LIB) -l$(HDF5LIBCPP) -lz

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECS = $(SOURCES:.cpp=)

all: $(EXECS)

$(EXECS): %: %.o
	$(CXX) $(CCOPTS) $< $(GTEST_MAINA) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CCOPTS) $< -c -o $@

gtest: $(EXECS)
	@(for exe in $(EXECS); do echo "Running $$exe"; ./$$exe --gtest_output=xml:$(ROOT_DIR)/$$exe.xml; done)

clean:
	@rm -f $(OBJECTS) $(EXECS)
