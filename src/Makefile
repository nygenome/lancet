BAMTOOLS_DIR := ../bamtools-2.3.0/
SEQLIB_DIR := ../SeqLib-1.1.1/

ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_SEQLIB_DIR := $(realpath $(SEQLIB_DIR))

CXX := g++
CXXFLAGS := -std=c++0x -Wno-deprecated -Wall -O3 -fexceptions
INCLUDES := -I$(ABS_SEQLIB_DIR) -I$(ABS_SEQLIB_DIR)/htslib -I$(ABS_BAMTOOLS_DIR)/include
LDFLAGS := -L$(ABS_SEQLIB_DIR)/bin -L$(ABS_BAMTOOLS_DIR)/lib
LDLIBS := -lz -ldl -lpthread -lcurl -lcrypto -lseqlib -lhts -lbamtools
 
all: lancet

lancet: Lancet.cc align.cc util.hh util.cc sha256.hh sha256.cc FET.hh ErrorCorrector.hh Mer.hh Ref.cc Ref.hh ReadInfo.hh ReadStart.hh Transcript.hh Variant.hh Variant.cc VariantDB.hh VariantDB.cc Edge.cc Edge.hh ContigLink.hh Node.cc Node.hh Path.cc Path.hh ContigLink.cc Graph.cc Graph.hh Microassembler.cc Microassembler.hh
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LDFLAGS) Lancet.cc Edge.cc Node.cc Graph.cc Microassembler.cc Ref.cc Path.cc ContigLink.cc align.cc util.cc sha256.cc VariantDB.cc Variant.cc -o lancet $(ABS_SEQLIB_DIR)/bin/libhts.a $(LDLIBS)

clean:
	rm -rf lancet
