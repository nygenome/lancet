BAMTOOLS_DIR := ../bamtools-2.3.0/
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
CXX := g++
CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
#CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -pg -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: Lancet

Lancet: Lancet.cc Lancet.hh align.cc util.cc Mer.hh Ref.cc Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh Variant.hh Variant.cc Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh
	$(CXX) $(CXXFLAGS) $(INCLUDES) Lancet.cc Edge.cc Node.cc Graph.cc Ref.cc Path.cc ContigLink.cc align.cc util.cc Variant.cc -o Lancet -lbamtools -lz