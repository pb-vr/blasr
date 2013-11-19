#
#
#
SHELL = /bin/bash
GCC = g++ -O3
# GCC = g++
SRC_DIR = src
BUILD_DIR = build
TARGET = pelusa

CCOPT = -g
LIBS = -lboost_program_options -lpthread -L ../../../ThirdParty/boost/boost_1_47_0/stage/lib/
# DEFINES = -DDB_SHORT -DREAD_36 -DMAXGAP=3 -DMAXHITS=10000
INCLUDES = -I ../../../ThirdParty/boost/boost_1_47_0/

OBJS = $(patsubst %.cpp,%.o, $(wildcard $(SRC_DIR)/*.cpp))

all: init $(TARGET)

init:
	@( [[ ! -e $(BUILD_DIR) ]] && mkdir $(BUILD_DIR) ) || true

$(TARGET): $(OBJS)
	$(GCC) -static $(DEFINES) $(INCLUDES) $(CCOPT) $(OBJS) -o $(BUILD_DIR)/$@ $(LIBS)

%.o: %.cpp
	$(GCC) $(DEFINES) $(INCLUDES) $(CCOPT) -c $< -o $@ $(LIBS)

INSTALL_DIR ?= .
install:
	/usr/bin/install -m 555 $(BUILD_DIR)/$(TARGET) $(INSTALL_DIR)/

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(OBJS)

tests:
	cram.py test/cram/*.t
