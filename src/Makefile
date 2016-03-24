BAMTOOLS_DIR := ../bamtools-2.3.0/
HTSLIB_DIR := ../htslib-1.3/

ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_HTSLIB_DIR := $(realpath $(HTSLIB_DIR))

CXX := g++
CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
#CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -pg -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -I$(ABS_HTSLIB_DIR)/ -L$(ABS_BAMTOOLS_DIR)/lib/ -L$(ABS_HTSLIB_DIR)/

all: Lancet

Lancet: Lancet.cc align.cc util.hh util.cc sha256.hh sha256.cc FET.hh ErrorCorrector.hh Mer.hh Ref.cc Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh CoverageWindow.hh Variant.hh Variant.cc VariantDB.hh VariantDB.cc Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh Microassembler.cc Microassembler.hh
	$(CXX) $(CXXFLAGS) $(INCLUDES) Lancet.cc Edge.cc Node.cc Graph.cc Microassembler.cc Ref.cc Path.cc ContigLink.cc align.cc util.cc sha256.cc VariantDB.cc Variant.cc -o Lancet $(ABS_HTSLIB_DIR)/libhts.a -lz -ldl -lbamtools -lpthread

clean:
	rm -rf Lancet;