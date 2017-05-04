#ifndef SEQLIB_FERMI_H
#define SEQLIB_FERMI_H

#include <string>
#include <cstdlib>
#include <iostream>

#include "SeqLib/BamRecord.h"

extern "C" 
{
#include "fermi-lite/htab.h"
#include "fermi-lite/fml.h"
#include "fermi-lite/bfc.h"
}

namespace SeqLib {

  /** Sequence assembly using FermiKit from Heng Li
   */
  class FermiAssembler {

  public:

    /** Create an empty FermiAssembler with default parameters */
    FermiAssembler ();

    /** Destroy by clearing all stored reads from memory */
    ~FermiAssembler();

    /** Provide a set of reads to be assembled 
     * @param brv Reads with or without quality scores
     * @note This will copy the reads and quality scores
     * into this object. Deallocation is automatic with object
     * destruction, or with ClearReads.
     */ 
    void AddReads(const BamRecordVector& brv);

    /** Clear all of the sequences and deallocate memory.
     * This is not required, as it will be done on object destruction
     */
    void ClearReads();

    /** Clear all of the contigs and deallocate memory.
     * This is not required, as it will be done on object destruction
     */
    void ClearContigs();

    /** Peform Bloom filter error correction of the reads
     * in place. */
    void CorrectReads();

    /** Peform Bloom filter error correction of the reads
     * in place. Also remove unique reads.
     */
    void CorrectAndFilterReads();
    
    /** Return the sequences in this object, which may have 
     * been error-corrected
     */
    UnalignedSequenceVector GetSequences() const;

    /** Perform the string graph assembly.
     * This will product the string graph,
     * and travserse the graph to emit contigs
     */
    void PerformAssembly();

    /** Return the assembled contigs
     * @return Assembled contigs in upper case strings (ACTGN)
     */
    std::vector<std::string> GetContigs() const;

    /** Perform assembly, without error correction */
    void DirectAssemble(float kcov);

    /** Set the minimum overlap between reads during string graph construction */
    void SetMinOverlap(uint32_t m) { opt.min_asm_ovlp = m; }

    /** Aggressively trim graph to discard heterozygotes. 
     * Suggested by lh3 for bacterial assembly
     * @note See: https://github.com/lh3/fermi-lite/blob/master/example.c
     */
    void SetAggressiveTrim() { opt.mag_opt.flag |= MAG_F_AGGRESSIVE; }

    /** From lh3: Drop an overlap if its length is below max_overlap * ratio
     * @param ratio Overlaps below ratio * max_overlap will be removed
     */
    void SetDropOverlapRatio(double ratio) { opt.mag_opt.min_dratio1 = ratio; }

    /** From lh3: Min k-mer & read count thresholds for ec and graph cleaning
     */
    void SetKmerMinThreshold(int min) { opt.min_cnt = min; }

    /** From lh3: Max k-mer & read count thresholds for ec and graph cleaning
     */
    void SetKmerMaxThreshold(int max) { opt.max_cnt = max; }

    // From lh3: retain a bubble if one side is longer than the other side by >INT-bp
    //void SetBubbleDifference(int bdiff) { opt.mag_opt.max_bdiff; }

    /** Return the minimum overlap parameter for this assembler */
    uint32_t GetMinOverlap() const { return opt.min_asm_ovlp; }

    /** Add a set of unaligned sequences to stage for assembly */
    void AddReads(const UnalignedSequenceVector& v);

    /** Add a single sequence to be assembled */
    void AddRead(const UnalignedSequence& r);

    /** Add a single sequence from an aligned reads to be assembled */
    void AddRead(const BamRecord& r);

    /** Return the number of sequences that are controlled by this assembler */
    size_t NumSequences() const { return n_seqs; }

  private:

    // reads to assemble
    fseq1_t *m_seqs;
  
    // size of m_seqs
    size_t m;
  
    std::vector<std::string> m_names;

    // number of base-pairs
    uint64_t size;
    
    // number of reads
    size_t n_seqs;

    // number of contigs
    int n_utg;
  
    // options
    fml_opt_t opt;

    // the unitigs
    fml_utg_t *m_utgs;

  };
  

}

#endif
