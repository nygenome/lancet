#ifndef SEQLIB_BFC_H
#define SEQLIB_BFC_H

extern "C" {
  #include "fermi-lite/bfc.h"
  #include "fermi-lite/fml.h"
}

#include "SeqLib/BamRecord.h"
#include "SeqLib/UnalignedSequence.h"

namespace SeqLib {

/** Class to perform error-correction using BFC algorithm
 *
 * BFC is designed and implemented by Heng Li (https://github.com/lh3/bfc). 
 * From Heng: It is a variant of the classical spectrum alignment algorithm introduced
 * by Pevzner et al (2001). It uses an exhaustive search to find a k-mer path 
 * through a read that minimizeds a heuristic objective function jointly considering
 * penalities on correction, quality and k-mer support.
 */
  class BFC {

  public:
    /** Construct a new BFC engine */
    BFC() {
      bfc_opt_init(&bfc_opt);
      ch = NULL;
      kmer = 0;
      flt_uniq = 0;
      n_seqs = 0;
      m_seqs = NULL;
      kcov = 0;
      tot_k = 0;
      sum_k = 0;
      tot_len = 0;
      m_seqs_size = 0;
    }

    ~BFC() {
      clear();
      if (ch)
	bfc_ch_destroy(ch);
    }

    /** Allocate a block of memory for the reads if the amount to enter is known 
     * @note This is not necessary, as reads will dynamically reallocate 
     */
    bool AllocateMemory(size_t n);

    /** Peform BFC error correction on the sequences stored in this object */
    bool ErrorCorrect();

    /** Train the error corrector using the reads stored in this object */
    bool Train();

    /** Add a sequence for either training or correction */
    bool AddSequence(const BamRecord& r);

    /** Add a sequence for either training or correction */
    bool AddSequence(const char* seq, const char* qual, const char* name);

    /** Set the k-mer size */
    void SetKmer(int k) { kmer = k; }

    /** Train error correction using sequences from aligned reads */
    void TrainCorrection(const BamRecordVector& brv);

    /** Train error correction from raw character strings */
    void TrainCorrection(const std::vector<char*>& v);

    /** Train and error correction on same reads */
    void TrainAndCorrect(const BamRecordVector& brv);

    /** Error correct a collection of reads */
    void ErrorCorrect(const BamRecordVector& brv);

    /** Error correct in place, modify sequence, and the clear memory from this object */
    void ErrorCorrectInPlace(BamRecordVector& brv);
    
    /** Error correct and add tag with the corrected sequence data, and the clear memory from this object 
     * @param brv Aligned reads to error correct
     * @param tag Tag to assign error corrected sequence to (eg KC)
     * @exception Throws an invalid_argument if tag is not length 2
     */
    void ErrorCorrectToTag(BamRecordVector& brv, const std::string& tag);
    
    /** Return the reads (error corrected if ran ErrorCorrect) */
    void GetSequences(UnalignedSequenceVector& v) const;

    /** Clear the stored reads */
    void clear();

    /** Filter reads with unique k-mers. Do after error correction */
    void FilterUnique();

    /** Return the calculated kcov */
    float GetKCov() const { return kcov; }

    /** Return the calculated kcov */
    int GetKMer() const { return kmer; }

    /** Return the number of sequences controlled by this */
    int NumSequences() const { return n_seqs; } 

  private:

    // the amount of memory allocated
    size_t m_seqs_size;

    void learn_correct();

    bfc_opt_t bfc_opt;

    // histogram of kmer occurences
    uint64_t hist[256];

    // diff histogram of kmers??
    uint64_t hist_high[64];

    uint64_t tot_len;

    uint64_t sum_k; // total valid kmer count (kmers above min_count) ? 

    // total number of kmers?
    uint64_t tot_k;

    //
    float kcov;

    // reads to correct in place
    fseq1_t * m_seqs;

    // number of sequeces
    size_t n_seqs;

    // fermi lite options
    fml_opt_t fml_opt;

    // vector of names
    std::vector<char*> m_names;

    // assign names, qualities and seq to m_seqs
    void allocate_sequences_from_reads(const BamRecordVector& brv);

    // assign names, qualities and seq to m_seqs
    void allocate_sequences_from_char(const std::vector<char*>& v);
    
    // do the actual read correction
    void correct_reads();

    // 0 turns off filter uniq
    int flt_uniq; // from fml_correct call
    
    int l_pre;

    // 0 is auto learn
    int kmer;

    // holds data after learning how to correct
    bfc_ch_t *ch;

    // holds data for actual error correction
    ec_step_t es;
  };

   }

#endif
