#ifndef SEQLIB_FASTQ_READER_H
#define SEQLIB_FASTQ_READER_H

#include <string>

#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>

// all kseq stuff is in UnaligedSequence
#include "SeqLib/UnalignedSequence.h"

namespace SeqLib{

  /** Simple reader for FASTA/FASTQ files */
class FastqReader {

 public:
  
  /** Construct an empty FASTQ/FASTA reader */
 FastqReader() {}
  
  /** Construct a reader and open a FASTQ/FASTA reader 
   * @param file Path to a FASTQ or FASTA file
   */
  FastqReader(const std::string& file);

  /** Open a FASTQ/FASTA file for reading
   * @param file Path to a FASTQ or FASTA file
   * @return Returns true if opening was successful
   */
  bool Open(const std::string& file);

  /** Retrieve the next sequence from the FASTA/FASTQ
   * @param s Sequence to be filled in with Name, Seq, Qual and Strand
   */
  bool GetNextSequence(UnalignedSequence& s);

  ~FastqReader() { 
    if (seq)  
      kseq_destroy(seq);
    if (fp)
      gzclose(fp);
  }


 private:
  
  std::string m_file;

  gzFile fp; // file handler for kseq
  kseq_t * seq; // current read

};

}

#endif
