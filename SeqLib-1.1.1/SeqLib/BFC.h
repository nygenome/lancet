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
      m_idx = 0;
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

    /** Peform BFC error correction on the sequences stored in this object */
    bool ErrorCorrect();

    /** Train the error corrector using the reads stored in this object */
    bool Train();

    /** Add a sequence for either training or correction 
     * @param seq A sequence to be copied into this object (A, T, C, G)
     */
    bool AddSequence(const char* seq, const char* qual, const char* name);

    /** Set the k-mer size for training 
     * @note zero is auto
     */
    void SetKmer(int k) { kmer = k; }

    /** Correct a single new sequence not stored in object 
     * @param str Sequence of string to correct (ACTG)
     * @param q Quality score of sequence to correct 
     * @value Returns true if corrected */
    bool CorrectSequence(std::string& str, const std::string& q);
    
    /** Clear the stored reads */
    void clear();

    /** Return the calculated kcov */
    float GetKCov() const { return kcov; }

    /** Return the calculated kcov */
    int GetKMer() const { return kmer; }

    /** Return the number of sequences controlled by this */
    int NumSequences() const { return n_seqs; } 

    /** Return the next sequence stored in object 
     * @param s Empty string to be filled.
     * @param q Empty string name to be filled.
     * @value True if string was filled with sequence. False if no more sequences.
     */
    bool GetSequence(std::string& s, std::string& q);

    /** Reset the sequence iterator inside GetSequence 
     */
    void ResetGetSequence() { m_idx = 0; };

  private:

    size_t m_idx;

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
    
    void allocate_sequences_from_strings(const std::vector<std::string>& v);
    
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
