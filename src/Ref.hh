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

	map<string,int> mertable;
	set<int> refcompids;

	int refnodes;
	int refcomp;
	int allcomp;

	bool indexed_m;
	
	vector<int> coverage; // k-mer coverage across the reference

	Ref_t(int k) : indexed_m(0) 
		{ K = k; }
	
	void setHdr(string hdr_) { hdr = hdr_; }	
	void setRawSeq(string rawseq_) { rawseq = rawseq_; }
	void setK(int k) { K = k; indexed_m = 0; mertable.clear(); resetCoverage(); }
	void setSeq(string seq_) { seq = seq_; coverage.resize(seq.size()); resetCoverage(); }

	void indexMers()
	{
		if (!indexed_m)
		{
			CanonicalMer_t cmer;
			for (unsigned int i = 0; i < seq.length() - K + 1; i++)
			{
				cmer.set(seq.substr(i, K));
			    mertable.insert(std::pair<string,int>(cmer.mer_m,0));
			}
			indexed_m = true;
		}
	}

	bool hasMer(const string & cmer)
	{
		indexMers();
		return mertable.count(cmer);
	}

	bool isRefComp(int comp)
	{
		return refcompids.find(comp) != refcompids.end();
	}
	
	void updateCoverage(const string & cmer) {
		indexMers();
		std::map<string,int>::iterator it = mertable.find(cmer);
		if (it != mertable.end()) {
			(*it).second += 1;
		}
	}
	
	// compute kmer coverage over the reference sequence 
	void computeCoverage() {
		CanonicalMer_t cmer;

		for (unsigned i = 0; i < seq.length() - K + 1; i++) 
		{
			cmer.set(seq.substr(i, K));			
			std::map<string,int>::iterator it = mertable.find(cmer.mer_m);
			if (it != mertable.end()) {
				int cov = (*it).second;
				if(i==0) {
					for (int j=i; j<K; j++) { coverage.at(j) = cov; }
				}
				else {
					coverage.at(i+K-1) = cov;
				}
			}		
		}
	}
	
	// print k-mer coverage along the reference
	void printKmerCoverage() {
	    cout << "cov: ";
		for (unsigned i=0; i<coverage.size(); i++) {
		    cout << " " << coverage.at(i);
		}
		cout << " len: " << coverage.size() << endl;
	}
	
	// reset coverage to 0
	void resetCoverage() {
		for (unsigned i=0; i<coverage.size(); i++) { 
			coverage.at(i) = 0;
		}
	}
	
	/*
	string getTrimSeq() {
		int ref_dist = trim3 - trim5 + K;
		return seq.substr(trim5, ref_dist);
	}
	*/
};

#endif
