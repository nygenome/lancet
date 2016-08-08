#ifndef REF_HH
#define REF_HH 1

/****************************************************************************
** Ref.hh
**
** Class for storing kmer information for the reference sequence
**
*****************************************************************************/

/************************** COPYRIGHT ***************************************
**
** New York Genome Center
**
** SOFTWARE COPYRIGHT NOTICE AGREEMENT
** This software and its documentation are copyright (2016) by the New York
** Genome Center. All rights are reserved. This software is supplied without
** any warranty or guaranteed support whatsoever. The New York Genome Center
** cannot be responsible for its use, misuse, or functionality.
**
** Version: 1.0.0
** Author: Giuseppe Narzisi
**
*************************** /COPYRIGHT **************************************/

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "Mer.hh"
#include "ReadInfo.hh"


using namespace std;

typedef struct cov_t
{
  int fwd; // total fwd coverage
  int rev; // total rev coverage 
  int minqv_fwd; // min base quality fwd coverage
  int minqv_rev; // min base quality rev coverage
  //int minmq_fwd; // min mapping quality fwd coverage
  //int minmq_rev; // min mapping quality rev coverage
} cov_t;

class Ref_t
{
public:

	unsigned short K;
	string hdr;
	string seq;

	string rawseq;

	string refchr;
	int refstart;
	int refend;

	unsigned short trim5;
	unsigned short trim3;

	// mapping of mers to fwd/rev counts (mer,cov_t)  
	//map<string,cov_t> mertable_nml;
	//map<string,cov_t> mertable_tmr;
	map<string,cov_t> * mertable_nml = NULL;
	map<string,cov_t> * mertable_tmr = NULL;
	
	set<int> refcompids;

	int refnodes;
	int refcomp;
	int allcomp;

	bool indexed_m;

	//vector<cov_t> normal_coverage; // normal k-mer coverage across the reference
	//vector<cov_t> tumor_coverage; // tumor k-mer coverage across the reference
	vector<cov_t> * normal_coverage = NULL; // normal k-mer coverage across the reference
	vector<cov_t> * tumor_coverage = NULL; // tumor k-mer coverage across the reference
	
	Ref_t(int k) : indexed_m(0) 
		{ K = k; }
	
	void setHdr(string hdr_) { hdr = hdr_; }
	void setRawSeq(string rawseq_) { rawseq = rawseq_; }
	void setK(int k) { K = k; indexed_m = 0; clear(); /*resetCoverage();*/ }
	void setSeq(string seq_) { seq = seq_; }
	//void setSeq(string seq_) { seq = seq_; normal_coverage.resize(seq.size()); tumor_coverage.resize(seq.size()); resetCoverage(); }

	void indexMers();
	bool hasMer(const string & cmer);
	bool isRefComp(int comp) { return refcompids.find(comp) != refcompids.end(); }
	
	void updateCoverage(const string & cmer, unsigned int strand, char sample);
	void computeCoverage(char sample);
	int getCovAt(unsigned pos, unsigned int strand, char sample);
	int getMinCovInKbp(unsigned pos, int K, char sample);
	void printKmerCoverage(char sample);
	void resetCoverage();
	void clear();
};

#endif
