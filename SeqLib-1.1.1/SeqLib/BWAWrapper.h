#ifndef SEQLIB_BWAWRAPPER_H
#define SEQLIB_BWAWRAPPER_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

#include "SeqLib/BamRecord.h"
#include "SeqLib/BamHeader.h"
#include "SeqLib/UnalignedSequence.h"

// all of the bwa and kseq stuff is in unaligned sequence
// best way I had to keep from clashes with klib macros

#define MEM_F_SOFTCLIP  0x200

namespace SeqLib {
 
/** Calls BWA-MEM on sequence queries and returns aligned reads, all in memory 
 * @note Calls core functions provided by Heng Li in BWA-MEM. https://github.com/lh3/bwa
 */
class BWAWrapper {

 public:

  /** Create an empty BWA MEM interface
   * @note Will initalize a BWA-MEM memopt structure
   * with the BWA-MEM defaults found in mem_opt_init.
   * Memory allocation and deallocation is automatically 
   * handled in constructor / destructor.
   */
  BWAWrapper() { 
    idx = 0;
    memopt = mem_opt_init();
    memopt->flag |= MEM_F_SOFTCLIP;
  }

  /** Destroy the BWAWrapper (deallocate index and options) */
  ~BWAWrapper() { 
    if (idx)
      bwa_idx_destroy(idx);
    if (memopt)
      free(memopt);
  }
  
  /** Retrieve the sequence name from its numeric ID 
   * @param id Numeric ID of the reference sequence
   * @exception throws an out_of_bounds if id not found
   */
  std::string ChrIDToName(int id) const;

  /** Create a BamHeader from the loaded index files */
  BamHeader HeaderFromIndex() const;

  /** Perform a BWA-MEM alignment of a single sequnece, and store hits in BamReadVector 
   * @param seq Sequence to be aligned
   * @param name Name of the sequence to be aligned
   * @param vec Alignment hits are appended to vec
   * @param hardclip Should the output BamRecord objects be hardclipped
   * @param keep_sec_with_frac_of_primary_score Set a threshold for whether a secondary alignment should be output
   * @param max_secondary Set a hard-limit on the number of secondary hits that will be reported
   */
  void AlignSequence(const std::string& seq, const std::string& name, BamRecordVector& vec, bool hardclip, 
			   double keep_sec_with_frac_of_primary_score, int max_secondary) const;

  /** Construct a new bwa index for this object. 
   * @param v vector of references to input (e.g. v = {{"r1", "AT"}};)
   * 
   * Throw an invalid_argument exception if any of the names or sequences
   * of the input UnalignedSequenceVector is empty
   */
  void ConstructIndex(const UnalignedSequenceVector& v);

  /** Retrieve a bwa index object from disk
   * @param file path a to an index fasta (index with bwa index)
   * @return True if successful
   * @note Will delete the old index if already stored
   */
  bool LoadIndex(const std::string& file);

  /** Dump the stored index to files
   * @note This does not write the fasta itself
   * @param index_name Write index files (*.sai, *.pac, *.ann, *.bwt, *.amb)
   * @return True if able to write index
   */
  bool WriteIndex(const std::string& index_name) const;

  /** Return the raw index in bwaidx_t form */
  bwaidx_t* GetIndex() const { return idx; }

  /** Return the number of reference sequences in current index
   * @return Number of reference sequences, or 0 if uninitialized
   */
  int NumSequences() const;

  /** Print some basic information about the loaded index */
  friend std::ostream& operator<<(std::ostream& out, const BWAWrapper& b);

  /** Set the gap open penalty
   * @param gap_open Gap open penalty. Default 6.
   * @exception Throws invalid_argument if gap_open < 0
   */
  void SetGapOpen(int gap_open);

  /** Set the gap open penalty
   * @param gap_ext Gap extension penalty. Default 1
   * @exception Throws invalid_argument if gap_ext < 0
   */
  void SetGapExtension(int gap_ext);

  /** Set the mismatch penalty
   * @param m Mismatch penalty (BWA-MEM b). Default 4
   * @exception Throws invalid_argument if m < 0
   */
  void SetMismatchPenalty(int m);

  /** Set the reseed trigger
   * @param r See BWA-MEM -r. Default 1.5
   * @exception Throws invalid_argument if r < 0
   */
  void SetReseedTrigger(float r);

  /** Set the SW alignment bandwidth
   * @param w See BWA-MEM -w. Default 100
   * @exception Throws invalid_argument if w < 0
   */
  void SetBandwidth(int w);

  /** Set the SW alignment Z dropoff
   * @param z See BWA-MEM -d. Default 100
   * @exception Throws invalid_argument if z < 0
   */
  void SetZDropoff(int z);

  /** Set the 3-prime clipping penalty
   * @param p See BWA-MEM -L. 
   * @exception Throws invalid_argument if p < 0
   */
  void Set3primeClippingPenalty(int p);

  /** Set the 5-prime clipping penalty
   * @param p See BWA-MEM -L. 
   * @exception Throws invalid_argument if p < 0
   */
  void Set5primeClippingPenalty(int p);

  /** Set the match score. Scales -TdBOELU
   * @note Since this scales penalty options, it should be 
   * probably be specified first, and then other options 
   * (eg gap penalty) can be Set explicitly afterwards.
   * @param a See BWA-MEM -A
   * @exception Throws invalid_argument if a < 0
   */
  void SetAScore(int a);

  /** Check if the index is empty */
  bool IsEmpty() const { return !idx; }
  
 private:

  // Construct a bam_hdr_t from a header string 
  bam_hdr_t* sam_hdr_read2(const std::string& hdr) const;

  // Store the options in memory
  mem_opt_t * memopt;

  // hold the full index structure
  bwaidx_t* idx;

  // Convert a bns to a header string 
  std::string bwa_print_sam_hdr2(const bntseq_t *bns, const char *hdr_line) const;

  // overwrite the bwa bwt_pac2pwt function
  bwt_t *seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr);

  // add an anns (chr annotation structure) 
  bntann1_t* seqlib_add_to_anns(const std::string& name, const std::string& seq, bntann1_t * ann, size_t offset);

  // overwrite the bwa-mem add1 function, which takes a sequence and adds to pac
  uint8_t* seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q);

  // make the pac structure (2-bit encoded packed sequence)
  uint8_t* seqlib_make_pac(const UnalignedSequenceVector& v, bool for_only);

  // write pac part of the index
  void seqlib_write_pac_to_file(const std::string& file) const;

  // write the bns file of the index
  std::string print_bns();
};

}


#endif
