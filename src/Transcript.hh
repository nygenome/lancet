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
	unsigned int ref_pos;
	char code;
	string ref;
	string qry;
	vector<int> cov_distr_N;
	vector<int> cov_distr_T;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative

	Transcript_t(int pos_, int ref_pos_, char code_, char ref_, char qry_, int covN_, int covT_, char prev_bp_ref_, char prev_bp_alt_)
		: pos(pos_), ref_pos(ref_pos_), code(code_)
	{ 
		ref = ref_;
		qry = qry_;
		cov_distr_N.push_back(covN_);
		cov_distr_T.push_back(covT_);
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	float getAvgCov(char sample) {
		float sum = 0;
		vector<int> cov_distr;
		
		if(sample == 'N') { cov_distr = cov_distr_N; }
		else if(sample == 'T') { cov_distr = cov_distr_T; }
		
		unsigned int n = cov_distr.size();
		for (unsigned int i = 0; i < n; i++)
		{
			sum += cov_distr[i];
		}
		return (float)sum/(float)n;
	}
	
	int getMinCov(char sample) {
		int min = 10000000;
		vector<int> cov_distr;
		
		if(sample == 'N') { cov_distr = cov_distr_N; }
		else if(sample == 'T') { cov_distr = cov_distr_T; }
		
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
