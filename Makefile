BAMTOOLS_DIR := ./bamtools-2.3.0/
HTSLIB_DIR := ./htslib-1.3/

ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
ABS_HTSLIB_DIR := $(realpath $(HTSLIB_DIR))

CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools htslib src

.PHONY : src
src:
	cd src; make; cd ../

.PHONY : bamtools
bamtools:
	mkdir $(ABS_BAMTOOLS_DIR)/build; cd $(ABS_BAMTOOLS_DIR)/build; cmake ..; make; cd ../../

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; make clean; cd ../../
	
.PHONY : htslib
bamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; cmake ..; make; cd ../../

#.PHONY : clean
clean:
	rm -rf $(ABS_BAMTOOLS_DIR)/build; rm -f src/Lancet; cd ..; cd ..;