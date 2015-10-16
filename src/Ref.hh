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

#include <string>
#include <set>

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

	set<string> mertable;
	set<int> refcompids;

	int refnodes;
	int refcomp;
	int allcomp;

	bool indexed_m;

	Ref_t(int k) : indexed_m(0) 
		{ K = k; }

	void setK(int k) { K = k; indexed_m = 0; mertable.clear(); }

	void indexMers()
	{
		if (!indexed_m)
		{
			CanonicalMer_t cmer;
			for (unsigned int i = 0; i < seq.length() - K + 1; i++)
			{
				cmer.set(seq.substr(i, K));
				mertable.insert(cmer.mer_m);
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
	
	/*
	string getTrimSeq() {
		int ref_dist = trim3 - trim5 + K;
		return seq.substr(trim5, ref_dist);
	}
	*/
};

#endif
