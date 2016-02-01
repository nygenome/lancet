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
	unsigned int end_pos;
	unsigned int ref_end_pos;
	string ref;
	string qry;
	bool isSomatic;

	float mean_cov_N;
	float mean_cov_T;
	float mean_cov_N_non0;
	float mean_cov_T_non0;
	
	float mean_ref_cov_N;
	float mean_ref_cov_T;
	float mean_ref_cov_N_non0;
	float mean_ref_cov_T_non0;
	
	int min_cov_N;
	int min_cov_T;
	int min_nonzero_cov_N;
	int min_nonzero_cov_T;
	
	int min_ref_cov_N;
	int min_ref_cov_T;
	int min_nonzero_ref_cov_N;
	int min_nonzero_ref_cov_T;
	
	vector<int> cov_distr_N;
	vector<int> cov_distr_T;
	vector<int> ref_cov_distr_N;
	vector<int> ref_cov_distr_T;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative

	Transcript_t(int pos_, int ref_pos_, char code_, char ref_, char qry_, int covN_, int covT_, int ref_covN_, int ref_covT_, char prev_bp_ref_, char prev_bp_alt_, int end_pos_, int ref_end_pos_, bool flag)
		: pos(pos_), ref_pos(ref_pos_), code(code_), end_pos(end_pos_), ref_end_pos(ref_end_pos_)
	{ 	
		isSomatic = flag;
		ref = ref_;
		qry = qry_;
		cov_distr_N.push_back(covN_);
		cov_distr_T.push_back(covT_);
		ref_cov_distr_N.push_back(ref_covN_);
		ref_cov_distr_T.push_back(ref_covT_);
		
		min_cov_N = covN_;
		min_cov_T = covT_;
		if(covN_ != 0) { min_nonzero_cov_N = ref_covN_; }
		if(covT_ != 0) { min_nonzero_cov_T = ref_covT_; }
		
		min_ref_cov_N = ref_covN_;
		min_ref_cov_T = ref_covT_;
		if(ref_covN_ != 0) { min_nonzero_ref_cov_N = ref_covN_; }
		if(ref_covT_ != 0) { min_nonzero_ref_cov_T = ref_covT_; }
		
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	// update coverage stats
	void updateStats() {
		float sumN = 0;
		float sumT = 0;
		float sumR = 0;
		float sumN_non0 = 0;
		float sumT_non0 = 0;
		float sumR_non0 = 0;
		unsigned int n;
		unsigned int n_non0;
			
		// update coverage stats for normal
		n = cov_distr_N.size();
		n_non0 = 0;
		for (unsigned int i = 0; i < n; i++) {
			sumN += cov_distr_N[i];
			if(cov_distr_N[i] != 0) { sumN_non0 += cov_distr_N[i]; n_non0++; }		
			if(cov_distr_N[i] < min_cov_N) { min_cov_N = cov_distr_N[i]; }
			if(cov_distr_N[i] < min_nonzero_cov_N && cov_distr_N[i] != 0) { min_nonzero_cov_N = cov_distr_N[i]; }
		}
		if(n>0) { mean_cov_N =(float)sumN/(float)n; }
		else { mean_cov_N = 0; }
		if (n_non0 > 0) { mean_cov_N_non0 =(float)sumN_non0/(float)n_non0; }
		else { mean_cov_N_non0 = 0; }
		
		// update coverage stats for tumor
		n = cov_distr_T.size();
		n_non0 = 0;
		for (unsigned int i = 0; i < n; i++) {
			sumT += cov_distr_T[i];
			if(cov_distr_T[i] != 0) { sumT_non0 += cov_distr_T[i]; n_non0++; }	
			if(cov_distr_T[i] < min_cov_T) { min_cov_T = cov_distr_T[i]; }
			if(cov_distr_T[i] < min_nonzero_cov_T && cov_distr_T[i] != 0) { min_nonzero_cov_T = cov_distr_T[i]; }
		}
		if(n>0) { mean_cov_T =(float)sumT/(float)n; }
		else { mean_cov_T = 0; }
		if (n_non0 > 0) { mean_cov_T_non0 =(float)sumT_non0/(float)n_non0; }
		else { mean_cov_T_non0 = 0; }
		
		// update coverage stats for reference (normal)
		n = ref_cov_distr_N.size();
		n_non0 = 0; sumR = 0;
		for (unsigned int i = 0; i < n; i++) {
			sumR += ref_cov_distr_N[i];
			if(ref_cov_distr_N[i] != 0) { sumR_non0 += ref_cov_distr_N[i]; n_non0++; }	
			if(ref_cov_distr_N[i] < min_ref_cov_N) { min_ref_cov_N = ref_cov_distr_N[i]; }
			if(ref_cov_distr_N[i] < min_nonzero_ref_cov_N && ref_cov_distr_N[i] != 0) { min_nonzero_ref_cov_N = ref_cov_distr_N[i]; }
		}
		if(n>0) { mean_ref_cov_N =(float)sumR/(float)n; }
		else { mean_ref_cov_N = 0; }
		if (n_non0 > 0) { mean_ref_cov_N_non0 =(float)sumR_non0/(float)n_non0; }
		else { mean_ref_cov_N_non0 = 0; }
		
		// update coverage stats for reference (normal)
		n = ref_cov_distr_T.size();
		n_non0 = 0; sumR = 0;
		for (unsigned int i = 0; i < n; i++) {
			sumR += ref_cov_distr_T[i];
			if(ref_cov_distr_T[i] != 0) { sumR_non0 += ref_cov_distr_T[i]; n_non0++; }	
			if(ref_cov_distr_T[i] < min_ref_cov_T) { min_ref_cov_T = ref_cov_distr_T[i]; }
			if(ref_cov_distr_T[i] < min_nonzero_ref_cov_T && ref_cov_distr_T[i] != 0) { min_nonzero_ref_cov_T = ref_cov_distr_T[i]; }
		}
		if(n>0) { mean_ref_cov_T =(float)sumR/(float)n; }
		else { mean_ref_cov_T = 0; }
		if (n_non0 > 0) { mean_ref_cov_T_non0 =(float)sumR_non0/(float)n_non0; }
		else { mean_ref_cov_T_non0 = 0; }
	}
	
	void addCovN(int c) { cov_distr_N.push_back(c); }	
	void addCovT(int c) { cov_distr_T.push_back(c); }
	void addRefCovN(int c) { ref_cov_distr_N.push_back(c); }	
	void addRefCovT(int c) { ref_cov_distr_T.push_back(c); }

	int getAvgCovN() { return ceil(mean_cov_N); }
	int getAvgCovT() { return ceil(mean_cov_T); }
	int getAvgNon0CovN() { return ceil(mean_cov_N_non0); }
	int getAvgNon0CovT() { return ceil(mean_cov_T_non0); }
	
	int getMinCovN() { return min_cov_N; }
	int getMinCovT() { return min_cov_T; }
	int getMinNon0CovN() { return min_nonzero_cov_N; }
	int getMinNon0CovT() { return min_nonzero_cov_T; }

	int getMinRefCovN() { return min_ref_cov_N; }
	int getMinRefCovT() { return min_ref_cov_T; }
	int getMinNon0RefCovN() { return min_nonzero_ref_cov_N; }
	int getMinNon0RefCovT() { return min_nonzero_ref_cov_T; }
	
	int getAvgRefCovN() { return ceil(mean_ref_cov_N); }
	int getAvgRefCovT() { return ceil(mean_ref_cov_T); }
	int getAvgNon0RefCovN() { return ceil(mean_ref_cov_N_non0); }
	int getAvgNon0RefCovT() { return ceil(mean_ref_cov_T_non0); }
};

#endif
