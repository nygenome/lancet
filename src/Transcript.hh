#ifndef TRANSCRIPT_HH
#define TRANSCRIPT_HH 1

/******************************************************************
** Transcript.hh
**
** Class for storing basic information about a genetic mutation 
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <string>

using namespace std;

class Transcript_t
{
public:

	unsigned int pos;
	char code;
	string ref;
	string qry;
	vector<int> cov_distr;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative

	Transcript_t(int pos_, char code_, char ref_, char qry_, int cov_, char prev_bp_ref_, char prev_bp_alt_)
		: pos(pos_), code(code_)
	{ 
		ref = ref_;
		qry = qry_;
		cov_distr.push_back(cov_);
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	float getAvgCov() {
		float sum = 0;
		unsigned int n = cov_distr.size();
		for (unsigned int i = 0; i < n; i++)
		{
			sum += cov_distr[i];
		}
		return (float)sum/(float)n;
	}
	
	int getMinCov() {
		int min = 10000000;
		unsigned int n = cov_distr.size();
		for (unsigned int i = 0; i < n; i++)
		{
			if (cov_distr[i] < min) { min = cov_distr[i]; }
		}
		if (min == 10000000) { min = -1; }
		return min;
	}

};

#endif
