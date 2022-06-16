BAMTOOLS_DIR := ./bamtools-2.5.2/
HTSLIB_DIR := ./htslib-1.15.1/

ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_HTSLIB_DIR := $(realpath $(HTSLIB_DIR))

CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools htslib lancet

.PHONY : lancet
lancet:
	cd src; make; cp lancet ../; cd ../

.PHONY : bamtools
bamtools:
	mkdir $(ABS_BAMTOOLS_DIR)/build; cd $(ABS_BAMTOOLS_DIR)/build; cmake -DCMAKE_INSTALL_PREFIX=../ ..; make; make install; cd ../../

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; make clean; cd ../../

.PHONY : htslib
htslib:
	cd $(ABS_HTSLIB_DIR); ./configure; make; cd ../

#.PHONY : clean
clean:
	 rm lancet src/lancet; rm -rf $(ABS_BAMTOOLS_DIR)/build; rm -rf $(ABS_BAMTOOLS_DIR)/include; rm -rf $(ABS_BAMTOOLS_DIR)/lib64; rm -rf $(ABS_BAMTOOLS_DIR)/bin; cd $(ABS_HTSLIB_DIR); make clean;
