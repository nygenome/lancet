/*
A significant portion of this code is derived from Heng Li's BWA
repository: https://github.com/lh3/bwa

BWA is copyrighted by Heng Li with the Apache2 License

*/

#include "SeqLib/BWAWrapper.h"

#include <stdexcept>
#include <sstream>
#include <iostream>

extern "C" {
  #include <string.h>
}

//#define DEBUG_BWATOOLS 1

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

namespace SeqLib {

  int BWAWrapper::NumSequences() const {
    
    if (!idx)
      return 0;

    return idx->bns->n_seqs;
    
  }

  std::string BWAWrapper::ChrIDToName(int id) const {

    if (!idx)
      throw std::runtime_error("Index has not be loaded / constructed");
    if (id < 0 || id >= idx->bns->n_seqs) 
      throw std::out_of_range("BWAWrapper::ChrIDToName - id out of bounds of refs in index for id of " + tostring(id) + " on IDX of size " + tostring(idx->bns->n_seqs));

    return std::string(idx->bns->anns[id].name);
  }

  BamHeader BWAWrapper::HeaderFromIndex() const 
  {

    std::string my_hdr = bwa_print_sam_hdr2(idx->bns, "");

    BamHeader hdr(my_hdr);
    //bam_hdr_t * hdr = bam_hdr_init();
    //bam_hdr_t * hdr = sam_hdr_read2(my_hdr); 
    //hdr->n_targets = idx->bns->n_seqs;
    //hdr->target_name = (char**)malloc(hdr->n_targets * sizeof(char*));
    //for (int i = 0; i < idx->bns->n_seqs; ++i) {
    //  hdr->target_name[i] = (char*)malloc( (strlen(idx->bns->anns[i].name) + 1) * sizeof(char));
    //  strcpy(hdr->target_name[i], idx->bns->anns[i].name);
    //}
    return hdr;
  }

  std::string BWAWrapper::bwa_print_sam_hdr2(const bntseq_t *bns, const char *hdr_line) const
  {
    std::string out;
    int i, n_SQ = 0;
    //extern char *bwa_pg;
    if (hdr_line) {
      const char *p = hdr_line;
      while ((p = strstr(p, "@SQ\t")) != 0) {
	if (p == hdr_line || *(p-1) == '\n') ++n_SQ;
	p += 4;
      }
    }

    //JEREMIAH
    // get makx size
    size_t max_s = 0;
    for (i = 0; i < bns->n_seqs; ++i)
      max_s = std::max(strlen(bns->anns[i].name), max_s);

    if (n_SQ == 0) {
      char buffer[max_s + 30];
      for (i = 0; i < bns->n_seqs; ++i) {
	//err_printf("@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	sprintf(buffer, "@SQ\tSN:%s\tLN:%d\n", bns->anns[i].name, bns->anns[i].len);
	out.append(buffer);
      }
    } else if (n_SQ != bns->n_seqs && bwa_verbose >= 2)
      fprintf(stderr, "[W::%s] %d @SQ lines provided with -H; %d sequences in the index. Continue anyway.\n", __func__, n_SQ, bns->n_seqs);
    
    if (hdr_line) { char buffer[200]; sprintf(buffer, "%s\n", hdr_line); out.append(buffer); } //err_printf("%s\n", hdr_line);
    //if (bwa_pg) { char buffer[100]; sprintf(buffer, "%s\n", bwa_pg); out.append(buffer); } // err_printf("%s\n", bwa_pg);
    
    return out;
  }
  
  
  
  void BWAWrapper::ConstructIndex(const UnalignedSequenceVector& v) {

    if (!v.size())
      return;

    // check the integrity of the input data
    for (UnalignedSequenceVector::const_iterator i = v.begin(); i != v.end(); ++i)
      if (i->Name.empty() || i->Seq.empty())
	throw std::invalid_argument("BWAWrapper::constructIndex - Reference sequences must have non-empty name and seq");
    
    if (idx) {
      std::cerr << "...clearing old index" << std::endl;
      bwa_idx_destroy(idx);
      idx = 0;
    }
    
    // allocate memory for idx
    idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;

    // construct the forward-only pac
    uint8_t* fwd_pac = seqlib_make_pac(v, true); //true->for_only

    // construct the forward-reverse pac ("packed" 2 bit sequence)
    uint8_t* pac = seqlib_make_pac(v, false); // don't write, becasue only used to make BWT

    size_t tlen = 0;
    for (UnalignedSequenceVector::const_iterator i = v.begin(); i != v.end(); ++i)
      tlen += i->Seq.length();
    
#ifdef DEBUG_BWATOOLS
    std::cerr << "ref seq length: " << tlen << std::endl;
#endif

    // make the bwt
    bwt_t *bwt;
    bwt = seqlib_bwt_pac2bwt(pac, tlen*2); // *2 for fwd and rev
    bwt_bwtupdate_core(bwt);
    free(pac); // done with fwd-rev pac 
    
    // construct sa from bwt and occ. adds it to bwt struct
    bwt_cal_sa(bwt, 32);
    bwt_gen_cnt_table(bwt);
        
    // make the bns
    bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
    bns->l_pac = tlen;
    bns->n_seqs = v.size();
    bns->seed = 11;
    bns->n_holes = 0;

    // make the anns
    bns->anns = (bntann1_t*)calloc(v.size(), sizeof(bntann1_t));
    size_t offset = 0;
    for (size_t k = 0; k < v.size(); ++k) {
      seqlib_add_to_anns(v[k].Name, v[k].Seq, &bns->anns[k], offset);
      offset += v[k].Seq.length();
    }

    //ambs is "holes", like N bases
    bns->ambs = 0; //(bntamb1_t*)calloc(1, sizeof(bntamb1_t));
    
    // make the in-memory idx struct
    idx->bwt = bwt;
    idx->bns = bns;
    idx->pac = fwd_pac;

    return;
    
  }

  void BWAWrapper::SetGapOpen(int gap_open) {
    if (gap_open < 0) 
      throw std::invalid_argument("BWAWrapper::SetGapOpen - gap_open must be >= zero");
    memopt->o_del = memopt->o_ins = gap_open;
  }

  void BWAWrapper::SetGapExtension(int gap_ext) {
    if (gap_ext < 0) 
      throw std::invalid_argument("BWAWrapper::SetGapExtension - gap extension must be >= zero");
    memopt->e_del = memopt->e_ins = gap_ext;
  }

  void BWAWrapper::SetMismatchPenalty(int m) {
    if (m < 0) 
      throw std::invalid_argument("BWAWrapper::SetMismatchPenalty - mismatch must be >= zero");
    memopt->b = m;

    bwa_fill_scmat(memopt->a, memopt->b, memopt->mat);
  }

  void BWAWrapper::SetZDropoff(int z) {
    if (z < 0) 
      throw std::invalid_argument("BWAWrapper::SetZDropoff - dropoff must be >= zero");
    memopt->zdrop = z;
  }

  void BWAWrapper::SetAScore(int a) {
    if (a < 0) 
      throw std::invalid_argument("BWAWrapper::SetAScore - dropoff must be >= zero");
    memopt->b *= a;
    memopt->T *= a;
    memopt->o_del *= a;
    memopt->o_ins *= a;
    memopt->e_del *= a;
    memopt->e_ins *= a;
    memopt->zdrop *= a;
    memopt->pen_clip5 *= a;
    memopt->pen_clip3 *= a;
    memopt->pen_unpaired *= a;
    memopt->a = a;
  }


  void BWAWrapper::Set3primeClippingPenalty(int p) {
    if (p < 0) 
      throw std::invalid_argument("BWAWrapper::Set3primeClippingPenalty - penalty must be >= zero");
    memopt->pen_clip3 = p;
  }

  void BWAWrapper::Set5primeClippingPenalty(int p) {
    if (p < 0) 
      throw std::invalid_argument("BWAWrapper::Set5primeClippingPenalty - penalty must be >= zero");
    memopt->pen_clip5 = p;
  }

  void BWAWrapper::SetBandwidth(int w) {
    if (w < 0) 
      throw std::invalid_argument("BWAWrapper::SetBandwidth - bandwidth must be >= zero");
    memopt->w = w;
  }

  void BWAWrapper::SetReseedTrigger(float r) {
    if (r < 0) 
      throw std::invalid_argument("BWAWrapper::SetReseedTrigger - reseed trigger must be >= zero");
    memopt->split_factor = r;
  }

  void BWAWrapper::AlignSequence(const std::string& seq, const std::string& name, BamRecordVector& vec, bool hardclip, 
				       double keep_sec_with_frac_of_primary_score, int max_secondary) const {
    
    // we haven't made an index, just return
    if (!idx)
      return;

    mem_alnreg_v ar;
    ar = mem_align1(memopt, idx->bwt, idx->bns, idx->pac, seq.length(), seq.data()); // get all the hits (was c_str())

#ifdef DEBUG_BWATOOLS
    std::cout << "num hits: " << ar.n << std::endl;
#endif    

    double primary_score = 0;

    int secondary_count = 0;

    //size_t num_secondary = 0;
    // loop through the hits
    for (size_t i = 0; i < ar.n; ++i) {

      if (ar.a[i].secondary >= 0 && (keep_sec_with_frac_of_primary_score < 0 || keep_sec_with_frac_of_primary_score > 1))
      	continue; // skip secondary alignments
      
      // get forward-strand position and CIGAR
      mem_aln_t a;

      a = mem_reg2aln(memopt, idx->bns, idx->pac, seq.length(), seq.c_str(), &ar.a[i]); 

      // if score not sufficient or past cap, continue
      bool sec_and_low_score =  ar.a[i].secondary >= 0 && (primary_score * keep_sec_with_frac_of_primary_score) > a.score;
      bool sec_and_cap_hit = ar.a[i].secondary >= 0 && (int)i > max_secondary;
      if (sec_and_low_score || sec_and_cap_hit) {
	free(a.cigar);
	continue;
      } else if (ar.a[i].secondary < 0) {
	primary_score = a.score;
	//num_secondary = 0;
      }

#ifdef DEBUG_BWATOOLS
      std::cerr << "allocing bamread" << std::endl;
#endif

      // instantiate the read
      BamRecord b;
      b.init();

      b.b->core.tid = a.rid;
      b.b->core.pos = a.pos;
      b.b->core.qual = a.mapq;
      b.b->core.flag = a.flag;
      b.b->core.n_cigar = a.n_cigar;
      
      // set dumy mate
      b.b->core.mtid = -1;
      b.b->core.mpos = -1;
      b.b->core.isize = 0;

      // if alignment is reverse, set it
      if (a.is_rev) 
	b.b->core.flag |= BAM_FREVERSE;

      std::string new_seq = seq;
      // if hardclip, figure out what to clip
      if (hardclip) {
	size_t tstart = 0;
	size_t len = 0;
	for (int i = 0; i < a.n_cigar; ++i) {
	  if (i == 0 && bam_cigar_op(a.cigar[i]) == BAM_CREF_SKIP) // first N (e.g. 20N50M)
	    tstart = bam_cigar_oplen(a.cigar[i]);
	  else if (bam_cigar_type(bam_cigar_op(a.cigar[i]))&1) // consumes query, but not N
	    len += bam_cigar_oplen(a.cigar[i]);
	}
	assert(len > 0);
	assert(tstart + len <= seq.length());
	new_seq = seq.substr(tstart, len);
      }

      // allocate all the data
      b.b->core.l_qname = name.length() + 1;
      b.b->core.l_qseq = new_seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
      b.b->l_data = b.b->core.l_qname + (a.n_cigar<<2) + ((b.b->core.l_qseq+1)>>1) + (b.b->core.l_qseq);
      b.b.get()->data = (uint8_t*)malloc(b.b.get()->l_data);

#ifdef DEBUG_BWATOOLS
      std::cerr << "memcpy" << std::endl;
#endif

      // allocate the qname
      memcpy(b.b->data, name.c_str(), name.length() + 1);

      // allocate the cigar. Reverse if aligned to neg strand, since mem_aln_t stores
      // cigars relative to referemce string oreiatnion, not forward alignment
      memcpy(b.b->data + b.b->core.l_qname, (uint8_t*)a.cigar, a.n_cigar<<2);

      // convert N to S or H
      int new_val = hardclip ? BAM_CHARD_CLIP : BAM_CSOFT_CLIP;
      uint32_t * cigr = bam_get_cigar(b.b);
      for (int k = 0; k < b.b->core.n_cigar; ++k) {
	if ( (cigr[k] & BAM_CIGAR_MASK) == BAM_CREF_SKIP) {
	  cigr[k] &= ~BAM_CIGAR_MASK;
	  cigr[k] |= new_val;
	}
      }
	
      // allocate the sequence
      uint8_t* m_bases = b.b->data + b.b->core.l_qname + (b.b->core.n_cigar<<2);

      // TODO move this out of bigger loop
      int slen = new_seq.length();
      int j = 0;
      if (a.is_rev) {
	for (int i = slen-1; i >= 0; --i) {
	  
	  // bad idea but works for now
	  // this is REV COMP things
	  uint8_t base = 15;
	  if (new_seq.at(i) == 'T')
	    base = 1;
	  else if (new_seq.at(i) == 'G')
	    base = 2;
	  else if (new_seq.at(i) == 'C')
	    base = 4;
	  else if (new_seq.at(i) == 'A')
	    base = 8;

	  m_bases[j >> 1] &= ~(0xF << ((~j & 1) << 2));   ///< zero out previous 4-bit base encoding
	  m_bases[j >> 1] |= base << ((~j & 1) << 2);  ///< insert new 4-bit base encoding
	  ++j;
	}
      } else {
	for (int i = 0; i < slen; ++i) {
	// bad idea but works for now
	  uint8_t base = 15;
	  if (new_seq.at(i) == 'A')
	    base = 1;
	  else if (new_seq.at(i) == 'C')
	    base = 2;
	  else if (new_seq.at(i) == 'G')
	    base = 4;
	  else if (new_seq.at(i) == 'T')
	    base = 8;
	  
	  m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
	  m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding

	}
      }

#ifdef DEBUG_BWATOOLS
      std::cerr << "memcpy3" << std::endl;
#endif

      // allocate the quality to NULL
      uint8_t* s = bam_get_qual(b.b);
      s[0] = 0xff;

      b.AddIntTag("NA", ar.n); // number of matches
      b.AddIntTag("NM", a.NM);

      if (a.XA)
	b.AddZTag("XA", std::string(a.XA));

      // add num sub opt
      b.AddIntTag("SB", ar.a[i].sub_n);
      b.AddIntTag("AS", a.score);

      // count num secondaries
      if (b.SecondaryFlag())
	++secondary_count;

      vec.push_back(b);

#ifdef DEBUG_BWATOOLS
      // print alignment
      printf("\t%c\t%s\t%ld\t%d\t", "+-"[a.is_rev], idx->bns->anns[a.rid].name, (long)a.pos, a.mapq);
      for (int k = 0; k < a.n_cigar; ++k) // print CIGAR
      	printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
      printf("\t%d\n", a.NM); // print edit distance
      std::cerr << "final done" << std::endl;
#endif
      
      free(a.cigar); // don't forget to deallocate CIGAR
    }
    free (ar.a); // dealloc the hit list

    // add the secondary counts
    for (BamRecordVector::iterator i = vec.begin(); i != vec.end(); ++i)
      i->AddIntTag("SQ", secondary_count);
    
}

  // modified from bwa (heng li)
uint8_t* BWAWrapper::seqlib_add1(const kseq_t *seq, bntseq_t *bns, uint8_t *pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q)
{
  bntann1_t *p;
  int lasts;
  if (bns->n_seqs == *m_seqs) {
    *m_seqs <<= 1;
    bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
  }
  p = bns->anns + bns->n_seqs;
  p->name = strdup((char*)seq->name.s);
  p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
  p->gi = 0; p->len = seq->seq.l;
  p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
  p->n_ambs = 0;
  for (size_t i = lasts = 0; i < seq->seq.l; ++i) {
    int c = nst_nt4_table[(int)seq->seq.s[i]];
    if (c >= 4) { // N
      if (lasts == seq->seq.s[i]) { // contiguous N
	++(*q)->len;
      } else {
	if (bns->n_holes == *m_holes) {
	  (*m_holes) <<= 1;
	  bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
	}
	*q = bns->ambs + bns->n_holes;
	(*q)->len = 1;
	(*q)->offset = p->offset + i;
	(*q)->amb = seq->seq.s[i];
	++p->n_ambs;
	++bns->n_holes;
      }
    }
    lasts = seq->seq.s[i];
    { // fill buffer
      if (c >= 4) c = lrand48()&3;
      if (bns->l_pac == *m_pac) { // double the pac size
	*m_pac <<= 1;
	pac = (uint8_t*)realloc(pac, *m_pac/4);
	memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
      }
      _set_pac(pac, bns->l_pac, c);
      ++bns->l_pac;
    }
  }
  ++bns->n_seqs;

  return pac;
}

  // modified from bwa (heng li)
uint8_t* BWAWrapper::seqlib_make_pac(const UnalignedSequenceVector& v, bool for_only)
{

  bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
  uint8_t *pac = 0; 
  int32_t m_seqs, m_holes;
  int64_t m_pac, l;
  bntamb1_t *q;

  bns->seed = 11; // fixed seed for random generator
  m_seqs = m_holes = 8; m_pac = 0x10000;
  bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
  bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
  pac = (uint8_t*) calloc(m_pac/4, 1);
  q = bns->ambs;

  // move through the unaligned sequences
  for (size_t k = 0; k < v.size(); ++k) {

    // make the ref name kstring
    kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
    name->l = v[k].Name.length() + 1;
    name->m = v[k].Name.length() + 3;
    name->s = (char*)calloc(name->m, sizeof(char));
    memcpy(name->s, v[k].Name.c_str(), v[k].Name.length()+1);
    
    // make the sequence kstring
    kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
    t->l = v[k].Seq.length();
    t->m = v[k].Seq.length() + 2;
    //t->s = (char*)calloc(v[k].Seq.length(), sizeof(char));
    t->s = (char*)malloc(t->m);
    memcpy(t->s, v[k].Seq.c_str(), v[k].Seq.length());
    
    // put into a kstring
    kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));  
    ks->seq = *t;
    ks->name = *name;
    
    // make the forward only pac
    pac = seqlib_add1(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);

    // clear it out
    free(name->s);
    free(name);
    free(t->s);
    free(t);
    //free(ks->name.s); 
    //free(ks->seq.s);
    //free(ks->f->buf);
    //free(
    free(ks);
    // NOTE free kstring_t?
    //kseq_destroy(s);
  }

  if (!for_only) 
    {
      // add the reverse complemented sequence
      m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
      pac = (uint8_t*)realloc(pac, m_pac/4);
      memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
      for (l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac)
	_set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
    }

  bns_destroy(bns);
  
  return pac;
}

  // modified from bwa (heng li)
bwt_t *BWAWrapper::seqlib_bwt_pac2bwt(const uint8_t *pac, int bwt_seq_lenr)
{

  bwt_t *bwt;
  ubyte_t *buf;
  int i;
  //FILE *fp;

  // initialization
  bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
  bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
  bwt->bwt_size = (bwt->seq_len + 15) >> 4;
  //fp = xopen(fn_pac, "rb");

  // prepare sequence
  //pac_size = (bwt->seq_len>>2) + ((bwt->seq_len&3) == 0? 0 : 1);
  //buf2 = (ubyte_t*)calloc(pac_size, 1);
  //err_fread_noeof(buf2, 1, pac_size, fp);
  //err_fclose(fp);
  memset(bwt->L2, 0, 5 * 4);
  buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
  for (i = 0; i < (int)bwt->seq_len; ++i) {
    buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
    ++bwt->L2[1+buf[i]];
  }
  for (i = 2; i <= 4; ++i) 
    bwt->L2[i] += bwt->L2[i-1];
  //free(buf2);

  // Burrows-Wheeler Transform
  bwt->primary = is_bwt(buf, bwt->seq_len);
  bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
  for (i = 0; i < (int)bwt->seq_len; ++i)
    bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
  free(buf);
  return bwt;
}

  // modified from bwa (heng li)
  bntann1_t* BWAWrapper::seqlib_add_to_anns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset) 
  {

    ann->offset = offset;
    ann->name = (char*)malloc(name.length()+1); // +1 for \0
    strncpy(ann->name, name.c_str(), name.length()+1);
    ann->anno = (char*)malloc(7);
    strcpy(ann->anno, "(null)\0");
    ann->len = seq.length();
    ann->n_ambs = 0; // number of "holes"
    ann->gi = 0; // gi?
    ann->is_alt = 0;
    
    return ann;
  }

  
  bool BWAWrapper::LoadIndex(const std::string& file)
  {

    // read in the bwa index
    bwaidx_t* idx_new = bwa_idx_load(file.c_str(), BWA_IDX_ALL);

    if (!idx_new) 
      return false;

    if (idx) {
      std::cerr << "...clearing old index" << std::endl;
      bwa_idx_destroy(idx);
    } 
    
    idx = idx_new;
    return true;
  }


  bool BWAWrapper::WriteIndex(const std::string& index_name) const
  {
    
    if (!idx) 
      return false;
    
    std::string bwt_name = index_name + ".bwt";
    std::string sa_name = index_name + ".sa";
    bwt_dump_bwt(bwt_name.c_str(), idx->bwt); 
    bwt_dump_sa(sa_name.c_str(), idx->bwt);
    bns_dump(idx->bns, index_name.c_str());
    seqlib_write_pac_to_file(index_name);

    return true;
  }

  // modified from bwa (heng li)
  void BWAWrapper::seqlib_write_pac_to_file(const std::string& file) const
  {
    // finalize .pac file
    FILE *fp;
    std::string nm = file + ".pac";
    fp = xopen(nm.c_str(), "wb");
    ubyte_t ct;
    err_fwrite(idx->pac, 1, (idx->bns->l_pac>>2) + ((idx->bns->l_pac&3) == 0? 0 : 1), fp);

    // the following codes make the pac file size always (l_pac/4+1+1)
    if (idx->bns->l_pac % 4 == 0) {
      ct = 0;
      err_fwrite(&ct, 1, 1, fp);
    }
    ct = idx->bns->l_pac % 4;
    err_fwrite(&ct, 1, 1, fp);

    // close .pac file
    err_fflush(fp);
    err_fclose(fp);
  }

  std::ostream& operator<<(std::ostream& out, const BWAWrapper& b) {
    out << "BNS: l_pac: " << b.idx->bns->l_pac << " n_seqs: " << b.idx->bns->n_seqs <<
      " seed: " << b.idx->bns->seed << " n_holes " << b.idx->bns->n_holes;
    return out;
  }
}

