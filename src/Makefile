BAMTOOLS_DIR := ../bamtools-2.5.1/
HTSLIB_DIR := ../htslib-1.8/

ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_HTSLIB_DIR := $(realpath $(HTSLIB_DIR))

0;95;cCXX := g++
CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib64/
#CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions -pg -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/ # compile for gprof
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/bamtools -I$(ABS_HTSLIB_DIR)/
LDFLAGS := -L$(ABS_BAMTOOLS_DIR)/lib64/ -L$(ABS_HTSLIB_DIR)/
LDLIBS := -llzma -lbz2 -lz -ldl -lpthread -lcurl -lcrypto -lbamtools

all: lancet

lancet: Lancet.cc Lancet.hh align.cc util.hh util.cc sha256.hh sha256.cc FET.hh ErrorCorrector.hh Mer.hh Ref.cc Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh Variant.hh Variant.cc VariantDB.hh VariantDB.cc Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh Microassembler.cc Microassembler.hh
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) Lancet.cc Edge.cc Node.cc Graph.cc Microassembler.cc Ref.cc Path.cc ContigLink.cc align.cc util.cc sha256.cc VariantDB.cc Variant.cc -o lancet $(ABS_HTSLIB_DIR)/libhts.a $(LDLIBS)

clean:
	rm -rf lancet;
