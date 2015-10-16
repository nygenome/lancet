BAMTOOLS_DIR := ./bamtools-2.3.0/
ABS_BAMTOOLS_DIR := $(realpath $(BAMTOOLS_DIR))
CXX := g++
CXXFLAGS := -Wno-deprecated -Wall -O3 -fexceptions -g -Wl,-rpath,$(ABS_BAMTOOLS_DIR)/lib/
INCLUDES := -I$(ABS_BAMTOOLS_DIR)/include/ -L$(ABS_BAMTOOLS_DIR)/lib/

all: bamtools src

.PHONY : src
src:
	cd src; make; cd ../

.PHONY : bamtools
bamtools:
	mkdir $(ABS_BAMTOOLS_DIR)/build; cd $(ABS_BAMTOOLS_DIR)/build; cmake ..; make; cd ../../

.PHONY : cleanbamtools
cleanbamtools:
	cd $(ABS_BAMTOOLS_DIR)/build; make clean; cd ../../

#.PHONY : clean
clean:
	rm -rf $(ABS_BAMTOOLS_DIR)/build; rm -f src/Lancet; cd ..; cd ..;