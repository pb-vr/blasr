.PHONY: all gtest clean

ROOT_DIR = ..
PBCPP_DIR= $(ROOT_DIR)/..
include $(ROOT_DIR)/common.mk

SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
EXECS = $(SOURCES:.cpp=)

SUB_DIRS = data datastructures files formatting \
		   algorithm utils

all: $(SOURCES) $(OBJECTS) $(EXECS)
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ]; \
		then $(MAKE) -C $$subdir -f Makefile; \
		fi;\
	done)

$(EXECS) : % : %.o
	$(CXX) $(CCOPTS) $< $(GTEST_MAINA) -o $@ $(LDFLAGS)

.cpp.o:
	$(CXX) $(CCOPTS) -c $< -o $@

gtest: $(EXECS)
	@(for exe in $(EXECS); do echo "Running $$exe"; ./$$exe; done;)
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ]; \
		then $(MAKE) -C $$subdir -f Makefile gtest; \
		fi;\
	done)

clean: 
	rm -f $(OBJECTS) $(EXECS)
	@(for subdir in $(SUB_DIRS); do \
		if [ -e $$subdir/Makefile ]; \
	    then $(MAKE) -C $$subdir -f Makefile clean; \
		fi;\
	done)

