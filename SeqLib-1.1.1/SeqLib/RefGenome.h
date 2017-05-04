#ifndef SEQLIB_REF_GENOME_H
#define SEQLIB_REF_GENOME_H

#include <string>
#include <cstdlib>
#include <iostream>

#include "htslib/htslib/faidx.h"

namespace SeqLib {
  
  /** Stores an indexed reference genome
   *
   * RefGenome is currently used as an interface to obtain
   * sequences from the reference given an interval.
   */
  class RefGenome {

  public:

    /** Create an empty RefGenome object */
    RefGenome() { index = NULL; }
    
    /** Destroy the malloc'ed faidx_t index inside object */
    ~RefGenome() { if (index) fai_destroy(index); }
    
    /** Query a region to get the sequence
     * @param chr_name name of the chr to query
     * @param p1 position 1. Zero-based
     * @param p2 position 2. Zero-based
     * 
     * @exception Throws an invalid_argument if p1 > p2, p1 < 0, p2 < 0, chr not found, or seq not found
     * @note This is currently NOT thread safe
     */
    std::string QueryRegion(const std::string& chr_name, int32_t p1, int32_t p2) const;

    /** Load an indexed reference sequence 
     * @param file Path to an indexed reference genome. See samtools faidx to create
     * @return True if succesfully loaded
     */
    bool LoadIndex(const std::string& file);
    
    /** Check if reference has been loaded */
    bool IsEmpty() const { 
      return (index == NULL); 
    }
    
  private:

    faidx_t * index;

  };
  

}

#endif
