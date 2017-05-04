#include "SeqLib/FastqReader.h"

#include <cctype>
#include <algorithm>

namespace SeqLib {

  bool FastqReader::Open(const std::string& f) {
    
    m_file = f;

    // check if file exists
    struct stat buffer;   
    if (stat(m_file.c_str(), &buffer) != 0) {
      std::cerr << "FastqReader: Failed to read non-existant file " << m_file << std::endl;
      return false;
    }

    fp = NULL;
    fp = (m_file != "-") ? gzopen(m_file.c_str(), "r") : gzdopen(fileno(stdin), "r");

    if (!fp) {
      std::cerr << "FastqReader: Failed to read " << m_file << std::endl;
      return false;
    }

    seq = kseq_init(fp); // set to first seq

    return true;
    
  }

  FastqReader::FastqReader(const std::string& file) : m_file(file) {
    Open(m_file);
  }

bool FastqReader::GetNextSequence(UnalignedSequence& s) {

  // kseq_read parses fastq and fasta

  if (!fp || !seq)
    return false;

  // no more reads
  if (kseq_read(seq) < 0)
    return false;

  if (seq->name.s)
    s.Name = std::string(seq->name.s, seq->name.l);
  if (seq->seq.s)
    s.Seq = std::string(seq->seq.s, seq->seq.l);
  if (seq->qual.s)
    s.Qual = std::string(seq->qual.s, seq->qual.l);
  
  return true;

}

}
