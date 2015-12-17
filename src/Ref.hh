#ifndef REF_HH
#define REF_HH 1

/******************************************************************
** Read.hh
**
** Class for storing kmer information for the reference
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <iostream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include "Mer.hh"


using namespace std;

class Ref_t
{
public:

	int K;
	string hdr;
	string seq;

	string rawseq;

	string refchr;
	int refstart;
	int refend;

	int trim5;
	int trim3;

	// mapping of mers to normal/tumor counts (mer,(n_cnt,t_cnt))  
	map<string,std::pair<int,int>> mertable;
	set<int> refcompids;

	int refnodes;
	int refcomp;
	int allcomp;

	bool indexed_m;
	
	vector<int> normal_coverage; // normal k-mer coverage across the reference
	vector<int> tumor_coverage; // tumor k-mer coverage across the reference

	Ref_t(int k) : indexed_m(0) 
		{ K = k; }
	
	void setHdr(string hdr_) { hdr = hdr_; }
	void setRawSeq(string rawseq_) { rawseq = rawseq_; }
	void setK(int k) { K = k; indexed_m = 0; mertable.clear(); resetCoverage(); }
	void setSeq(string seq_) { seq = seq_; normal_coverage.resize(seq.size()); tumor_coverage.resize(seq.size()); resetCoverage(); }

	void indexMers();
	bool hasMer(const string & cmer);
	bool isRefComp(int comp) { return refcompids.find(comp) != refcompids.end(); }
	
	void updateCoverage(const string & cmer, char sample);
	void computeCoverage();
	int getCovAt(unsigned pos, char sample);
	void printKmerCoverage(char sample);
	void resetCoverage();
};

#endif
