/*
A significant portion of this code is derived from Heng Li's BFC
repository: https://github.com/lh3/bfc

BFC is copyrighted by Heng Li with the following license:

The MIT License
 
Copyright (c) 2015 Broad Institute
 
Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:
 
The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.
 
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "SeqLib/BFC.h"

#include <stdexcept>
#include <algorithm>

namespace SeqLib {

  bool BFC::AllocateMemory(size_t n) {

    if (n <= 0)
      return false;
    
    m_seqs_size = n;
    m_seqs = (fseq1_t*)realloc(m_seqs, n * sizeof(fseq1_t));

    if (!m_seqs)
      return false;

    return true;
  }

  bool BFC::AddSequence(const BamRecord& r) {

    //char* s = strdup(r.Sequence().c_str());
    const char* q = bam_get_qname(r.raw());
    //uint8_t* l = bam_get_qual(r.raw());
    //char* qual = (char*)malloc(r.Length() + 1);
    //if (l)
    //  for (size_t i = 0; i < r.Length(); ++i)
    //	qual[i] = l[i] + 33;
    //qual[r.Length()] = '\0';
    
    bool ret = AddSequence(r.Sequence().c_str(), r.Qualities().c_str(), q);

    //if (s)
    //  free(s);
    //if (qual)
    //  free(qual);

    return ret;

  }

  bool BFC::AddSequence(const char* seq, const char* qual, const char* name) {

    // do the intial allocation
    if (n_seqs == 0 && !m_seqs) {
      m_seqs_size = 32;
      m_seqs = (fseq1_t*)malloc(m_seqs_size * sizeof(fseq1_t));
    }
    // realloc if not enough space
    else if (n_seqs >= m_seqs_size) {
      m_seqs_size = 2 * m_seqs_size;
      m_seqs = (fseq1_t*)realloc(m_seqs, m_seqs_size * sizeof(fseq1_t));
    }

    if (!m_seqs)
      return false;

    // make sure seq and qual are even valid (if qual provided)
    if (strlen(qual) && seq && qual) 
      if (strlen(seq) != strlen(qual))
	return false;
    if (!strlen(seq))
      return false;

    fseq1_t *s;
    
    s = &m_seqs[n_seqs];
    
    s->seq   = strdup(seq);
    s->qual = 0;
    if (strlen(qual)) {
      s->qual  = strdup(qual);
    }
    
    s->l_seq = strlen(seq);
    n_seqs++;
    
    m_names.push_back(strdup(name));

    assert(m_names.size() == n_seqs);

    return true;
  }


  bool BFC::ErrorCorrect() {
    correct_reads();
    return true;
  }

  bool BFC::Train() {
    learn_correct();
    return true;
  }

  void BFC::TrainAndCorrect(const BamRecordVector& brv) {

    // if already allocated, clear the old ones
    clear();

    // send reads to string
    allocate_sequences_from_reads(brv);

    // learn how to correct
    learn_correct();

    // do the correction
    correct_reads();

  }

  void BFC::TrainCorrection(const std::vector<char*>& v) {

    // if already allocated, clear the old ones
    clear();

    // set m_seqs and n_seqs
    allocate_sequences_from_char(v);

    // learn correct, set ch
    learn_correct();


  }
  
  void BFC::TrainCorrection(const BamRecordVector& brv) {

    // if already allocated, clear the old ones
    clear();

    // set m_seqs and n_seqs
    allocate_sequences_from_reads(brv);

    // learn correct, set ch
    learn_correct();
  }

  void BFC::ErrorCorrectToTag(BamRecordVector& brv, const std::string& tag) {
    
    if (tag.length() != 2)
      throw std::invalid_argument("Tag length should be 2");

    flt_uniq = 0;

    // if already allocated, clear the old ones
    clear();

    // send reads to string
    allocate_sequences_from_reads(brv);

    // do the correction
    correct_reads();

    assert(n_seqs == brv.size());
    for (size_t i = 0; i < n_seqs; ++i) {
      std::string str = std::string(m_seqs[i].seq);
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      brv[i].AddZTag("KC", str);
    }
    
    clear();

  }

  void BFC::ErrorCorrect(const BamRecordVector& brv) {

    flt_uniq = 0;

    // if already allocated, clear the old ones
    clear();

    // send reads to string
    allocate_sequences_from_reads(brv);

    // do the correction
    correct_reads();
  }

  void BFC::ErrorCorrectInPlace(BamRecordVector& brv) {

    flt_uniq = 0;
    clear();

    allocate_sequences_from_reads(brv);

    correct_reads();
    
    assert(n_seqs == brv.size());
    for (size_t i = 0; i < n_seqs; ++i) {
      std::string str = std::string(m_seqs[i].seq);
      std::transform(str.begin(), str.end(),str.begin(), ::toupper);
      brv[i].SetSequence(str);
    }

    clear();
  }
    
  void BFC::GetSequences(UnalignedSequenceVector& v) const {

    for (size_t i = 0; i < n_seqs; ++i)
      if (m_seqs[i].seq) { // wont be here if filter unique was called
	std::string str = std::string(m_seqs[i].seq);
	std::transform(str.begin(), str.end(),str.begin(), ::toupper);
	std::string name = m_names[i] ? std::string(m_names[i]) : std::string();
	std::string qual = m_seqs[i].qual ? std::string(m_seqs[i].qual) : std::string();
	v.push_back(UnalignedSequence(name, str, qual));
      }
    
  }

  void BFC::allocate_sequences_from_char(const std::vector<char*>& v) {

    m_seqs_size = v.size();
    m_seqs = (fseq1_t*)malloc(v.size() * sizeof(fseq1_t));
    
    uint64_t size = 0;
    for (std::vector<char*>::const_iterator r = v.begin(); r != v.end(); ++r) {
    //    for (auto& r : v) {
      fseq1_t *s;
      
      s = &m_seqs[n_seqs];
      
      s->seq   = strdup(*r);
      s->qual  = NULL; 
      
      s->l_seq = strlen(*r);
      size += m_seqs[n_seqs++].l_seq;
    }
    return;

  }

  void BFC::allocate_sequences_from_reads(const BamRecordVector& brv) {
      
    // alloc the memory
    m_seqs_size = brv.size();
    m_seqs = (fseq1_t*)malloc(brv.size() * sizeof(fseq1_t));
    
    uint64_t size = 0;
    for (BamRecordVector::const_iterator r = brv.begin(); r != brv.end(); ++r) {
      //    for (auto& r : brv) {
      m_names.push_back(strdup(r->Qname().c_str()));

      fseq1_t *s;
      
      s = &m_seqs[n_seqs];

      std::string qs = r->QualitySequence();
      s->seq   = strdup(qs.c_str());
      s->qual  = strdup(r->Qualities().c_str());
      
      s->l_seq = qs.length();
      size += m_seqs[n_seqs++].l_seq;
    }
    return;
  }

  void free_char(char*& c) {
    if (c) {
      free (c);
      c = NULL;
    }
  }

  void BFC::clear() {
    
    assert(m_names.size() == n_seqs);
    for (size_t i = 0; i < n_seqs; ++i) {
      free_char(m_names[i]);
      free_char(m_seqs[i].seq);
      free_char(m_seqs[i].qual);
    }

    if (m_seqs)
      free(m_seqs);
    m_seqs = 0;
    n_seqs = 0;

    m_names.clear();
    m_seqs_size = 0;

  }

  void BFC::learn_correct() {
    
    // options
    fml_opt_init(&fml_opt);
    
    // if kmer is 0, fix 
    if (kmer <= 0) {
      fml_opt_adjust(&fml_opt, n_seqs, m_seqs);
      kmer = fml_opt.ec_k;
    }

    // initialize BFC options
    for (size_t i = 0; i < n_seqs; ++i) 
      tot_len += m_seqs[i].l_seq; // compute total length
    bfc_opt.l_pre = tot_len - 8 < 20? tot_len - 8 : 20;
    
    //  setup the counting of kmers
    memset(&es, 0, sizeof(ec_step_t));
    //kmer is learned before this
    
    bfc_opt.k = kmer;
    
    //es.opt = &bfc_opt, es.n_seqs = n_seqs, es.seqs = m_seqs, es.flt_uniq = flt_uniq;
    
    // hold count info. also called bfc_ch_s. Composed of
    //    int k
    //    int l_pre
    //    cnthash_t **h
    //        h is of size 1<<l_pre (2^l_pre). It is array of hash tables
    //        h[i] is initialized with kh_init(cnt) which makes a cnthash_t
    // bfc_ch_t *ch; // set in BFC.h
    
    // do the counting
    ch = fml_count(n_seqs, m_seqs, bfc_opt.k, bfc_opt.q, bfc_opt.l_pre, bfc_opt.n_threads);

#ifdef DEBUG_BFC
    // size of random hash value
    khint_t k;
    int* ksize = (int*)calloc(1<<ch->l_pre, sizeof(int));
    for (int i = 0; i < (1<<ch->l_pre); ++i) {
      for (k = kh_begin(ch->h[i]); k != kh_end(ch->h[i]); ++k)
        ++ksize[i];
      fprintf(stderr, "K: %d S: %d\n", i, ksize[i]);
    }
#endif
  }

  void BFC::correct_reads() {
    
    assert(kmer > 0);

    es.ch = ch;
    es.opt = &bfc_opt;
    es.n_seqs = n_seqs;
    es.seqs = m_seqs;
    es.flt_uniq = flt_uniq;

    // make the histogram?
    // es.ch is unchanged (const)
    int mode = bfc_ch_hist(es.ch, hist, hist_high);

    for (int i = fml_opt.min_cnt; i < 256; ++i) 
      sum_k += hist[i], tot_k += i * hist[i];    

#ifdef DEBUG_BFC
    std::cerr << " sum_k " << sum_k << " tot_k " << tot_k << std::endl;
    fprintf(stderr, "MODE: %d\n", mode);
    for (int i = fml_opt.min_cnt; i < 256; ++i) {
      fprintf(stderr, "hist[%d]: %d\n",i,hist[i]);
    }
    for (int i = fml_opt.min_cnt; i < 64; ++i) {
      fprintf(stderr, "hist_high[%d]: %d\n",i,hist_high[i]);
    }
#endif
    
    kcov = (float)tot_k / sum_k;
    bfc_opt.min_cov = (int)(BFC_EC_MIN_COV_COEF * kcov + .499);
    bfc_opt.min_cov = bfc_opt.min_cov < fml_opt.max_cnt? bfc_opt.min_cov : fml_opt.max_cnt;
    bfc_opt.min_cov = bfc_opt.min_cov > fml_opt.min_cnt? bfc_opt.min_cov : fml_opt.min_cnt;

#ifdef DEBUG_BFC
    fprintf(stderr, "kcov: %f mincov: %d  mode %d \n", kcov, bfc_opt.min_cov, mode);  
#endif
  
    // do the actual error correction
    kmer_correct(&es, mode, ch);

    return;


  }

    void BFC::FilterUnique() {
      flt_uniq = 1;
      correct_reads();

      size_t count = 0;
      for (size_t i = 0; i < n_seqs; ++i)
	if (m_seqs[i].seq)
	  ++count;
    }
  
}
