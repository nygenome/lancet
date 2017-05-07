##INCLUDES=-I/xchip/gistic/Jeremiah/software/seqan-library-2.0.2/include -I /xchip/gistic/Jeremiah/GIT/SeqLib/src -I/xchip/gistic/Jeremiah/GIT/SeqLib/htslib -I/xchip/gistic/Jeremiah/software/boost_1.61.0_gcc5.1 -I/xchip/gistic/Jeremiah/software/bamtools-2.4.0/include
INCLUDES=-I/xchip/gistic/Jeremiah/software/seqan-library-2.2.0/include -I /xchip/gistic/Jeremiah/GIT/SeqLib -I/xchip/gistic/Jeremiah/GIT/SeqLib/htslib -I/xchip/gistic/Jeremiah/software/boost_1.61.0_gcc5.1 -I/xchip/gistic/Jeremiah/software/bamtools-2.4.0/include -I/xchip/gistic/Jeremiah/software/bzip2-1.0.6
LIBS=/xchip/gistic/Jeremiah/software/bamtools-2.4.0/lib/libbamtools.a /xchip/gistic/Jeremiah/GIT/SeqLib/src/libseqlib.a /xchip/gistic/Jeremiah/GIT/SeqLib/htslib/libhts.a /xchip/gistic/Jeremiah/software/boost_1.61.0_gcc5.1/stage/lib/libboost_timer.a /xchip/gistic/Jeremiah/software/boost_1.61.0_gcc5.1/stage/lib/libboost_chrono.a /xchip/gistic/Jeremiah/software/boost_1.61.0_gcc5.1/stage/lib/libboost_system.a /xchip/gistic/Jeremiah/software/bzip2-1.0.6/libbz2.a
CXXFLAGS=-W -Wall -pedantic -std=c++14 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0

binaries=benchmark

all: benchmark.o
	g++ benchmark.o -o benchmark $(LIBS) -lrt -lpthread -lz -lm

benchmark.o: benchmark.cpp
	g++ -c benchmark.cpp $(INCLUDES) $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(binaries) *.o
