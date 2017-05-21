BAMTOOLS_DIR := ./bamtools-2.3.0/
SEQLIB_DIR := ./SeqLib-1.1.1/

ABS_LANCET_DIR := $(realpath .)
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_SEQLIB_DIR := $(realpath $(SEQLIB_DIR))

CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools seqlib lancet

.PHONY : lancet
lancet:
	cd src && make && cp lancet $(ABS_LANCET_DIR) && cd $(ABS_LANCET_DIR)

.PHONY : bamtools
bamtools:
	cd $(ABS_BAMTOOLS_DIR) && mkdir -p build && cd build && cmake .. && make && cd $(ABS_LANCET_DIR)

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR) && rm -rf build include lib bin && cd $(ABS_LANCET_DIR)

.PHONY : seqlib
seqlib:
	cd $(ABS_SEQLIB_DIR)/htslib && ./configure --enable-libcurl && make
	cd $(ABS_SEQLIB_DIR) && ./configure LDFLAGS="-lcurl -lcrypto" && make && make install && cd $(ABS_LANCET_DIR)

.PHONY : cleanseqlib
cleanseqlib:
	cd $(ABS_SEQLIB_DIR)/htslib && make clean
	cd $(ABS_SEQLIB_DIR)/bwa && make clean
	cd $(ABS_SEQLIB_DIR)/fermi-lite && make clean
	find $(ABS_SEQLIB_DIR) -name Makefile | grep -v "bwa\|benchmark\|fermi\|htslib" | xargs rm -f
	cd $(ABS_SEQLIB_DIR) && rm -rf bin && cd $(ABS_LANCET_DIR)

.PHONY : clean
clean:
	rm -f lancet src/lancet
	make cleanbamtools
	make cleanseqlib
