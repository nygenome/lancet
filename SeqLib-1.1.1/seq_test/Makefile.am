bin_PROGRAMS = seq_test

seq_test_CPPFLAGS = \
     -I../src \
     -I../htslib \
     -I.. \
     -I../fermi-lite --coverage

seq_test_LDADD = \
	../fermi-lite/libfml.a \
	../bwa/libbwa.a \
	../htslib/libhts.a \
	-lboost_unit_test_framework -lboost_system -lboost_timer -lboost_chrono

seq_test_LDFLAGS = --coverage -pthread

seq_test_SOURCES = seq_test.cpp \
	../src/BFC.cpp ../src/GenomicRegion.cpp \
	../src/BamWriter.cpp ../src/BamReader.cpp \
	../src/ReadFilter.cpp ../src/BamRecord.cpp \
	../src/BWAWrapper.cpp \
        ../src/RefGenome.cpp ../src/SeqPlot.cpp ../src/BamHeader.cpp \
	../src/FermiAssembler.cpp ../src/ssw_cpp.cpp ../src/ssw.c ../src/jsoncpp.cpp
