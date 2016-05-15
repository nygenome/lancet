#ifndef TRANSCRIPT_HH
#define TRANSCRIPT_HH 1

/****************************************************************************
** Transcript.hh
**
** Class for storing basic information about a genetic mutation 
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

#include <string>
#include <algorithm>

using namespace std;

class Transcript_t
{
public:

	unsigned int pos;
	unsigned int ref_pos;
	unsigned int start_pos;
	char code;
	unsigned int end_pos;
	unsigned int ref_end_pos;
	string ref;
	string qry;
	bool isSomatic;

	float mean_cov_N_fwd;
	float mean_cov_N_rev;
	float mean_cov_T_fwd;
	float mean_cov_T_rev;
	
	float mean_cov_N_non0_fwd;
	float mean_cov_N_non0_rev;
	float mean_cov_T_non0_fwd;
	float mean_cov_T_non0_rev;
	
	float mean_ref_cov_N;
	float mean_ref_cov_T;
	float mean_ref_cov_N_non0;
	float mean_ref_cov_T_non0;
	
	int min_cov_N_fwd;
	int min_cov_N_rev;
	int min_cov_T_fwd;
	int min_cov_T_rev;
	
	int median_cov_N_fwd;
	int median_cov_N_rev;
	int median_cov_T_fwd;
	int median_cov_T_rev;
	
	int min_nonzero_cov_N_fwd;
	int min_nonzero_cov_N_rev;
	int min_nonzero_cov_T_fwd;
	int min_nonzero_cov_T_rev;
	
	int min_ref_cov_N;
	int min_ref_cov_T;
	int median_ref_cov_N;
	int median_ref_cov_T;
	int min_nonzero_ref_cov_N;
	int min_nonzero_ref_cov_T;
	
	vector<int> cov_distr_N_fwd;
	vector<int> cov_distr_T_fwd;
	vector<int> cov_distr_N_rev;
	vector<int> cov_distr_T_rev;
	vector<int> ref_cov_distr_N;
	vector<int> ref_cov_distr_T;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative

	Transcript_t(int pos_, int ref_pos_, int start_pos_, char code_, char ref_, char qry_, int covN_fwd_, int covN_rev_, int covT_fwd_, int covT_rev_, int ref_covN_, int ref_covT_, char prev_bp_ref_, char prev_bp_alt_, int end_pos_, int ref_end_pos_, bool flag)
		: pos(pos_), ref_pos(ref_pos_), start_pos(start_pos_), code(code_), end_pos(end_pos_), ref_end_pos(ref_end_pos_)
	{ 	
		isSomatic = flag;
		ref = ref_;
		qry = qry_;
		cov_distr_N_fwd.push_back(covN_fwd_);
		cov_distr_N_rev.push_back(covN_rev_);
		cov_distr_T_fwd.push_back(covT_fwd_);
		cov_distr_T_rev.push_back(covT_rev_);
		ref_cov_distr_N.push_back(ref_covN_);
		ref_cov_distr_T.push_back(ref_covT_);
		
		min_cov_N_fwd = covN_fwd_;
		min_cov_N_rev = covN_rev_;
		min_cov_T_fwd = covT_fwd_;
		min_cov_T_rev = covT_rev_;

		min_nonzero_cov_N_fwd = 0;
		min_nonzero_cov_N_rev = 0;
		min_nonzero_cov_T_fwd = 0;
		min_nonzero_cov_T_rev = 0;
		if(min_cov_N_fwd != 0) { min_nonzero_cov_N_fwd = min_cov_N_fwd; }
		if(min_cov_N_rev != 0) { min_nonzero_cov_N_rev = min_cov_N_rev; }
		if(min_cov_T_fwd != 0) { min_nonzero_cov_T_fwd = min_cov_T_fwd; }
		if(min_cov_T_rev != 0) { min_nonzero_cov_T_rev = min_cov_T_rev; }
				
		min_ref_cov_N = ref_covN_;
		min_ref_cov_T = ref_covT_;

		min_nonzero_ref_cov_N = 0;
		min_nonzero_ref_cov_T = 0;				
		if(ref_covN_ != 0) { min_nonzero_ref_cov_N = ref_covN_; }
		if(ref_covT_ != 0) { min_nonzero_ref_cov_T = ref_covT_; }
		
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	// update coverage stats for normal, tumor and the reference
	void updateStats() {				

		// normal
		computeStats(cov_distr_N_fwd, min_cov_N_fwd, min_nonzero_cov_N_fwd, mean_cov_N_fwd, mean_cov_N_non0_fwd, median_cov_N_fwd);
		computeStats(cov_distr_N_rev, min_cov_N_rev, min_nonzero_cov_N_rev, mean_cov_N_rev, mean_cov_N_non0_rev, median_cov_N_rev);
		
		//tumor
		computeStats(cov_distr_T_fwd, min_cov_T_fwd, min_nonzero_cov_T_fwd, mean_cov_T_fwd, mean_cov_T_non0_fwd, median_cov_T_fwd);
		computeStats(cov_distr_T_rev, min_cov_T_rev, min_nonzero_cov_T_rev, mean_cov_T_rev, mean_cov_T_non0_rev, median_cov_T_rev);
		
		// normal and tumor reference coverage
		computeStats(ref_cov_distr_N, min_ref_cov_N, min_nonzero_ref_cov_N, mean_ref_cov_N, mean_ref_cov_N_non0, median_ref_cov_N);
		computeStats(ref_cov_distr_T, min_ref_cov_T, min_nonzero_ref_cov_T, mean_ref_cov_T, mean_ref_cov_T_non0, median_ref_cov_T);
	}
	
	// compute coverage stats
	void computeStats(vector<int> &cov_distr, int &min_cov, int &min_nonzero_cov, float &mean_cov, float &mean_cov_non0, int &median_cov) {
		float sum= 0;
		unsigned int n = cov_distr.size();
		unsigned int n_non0 = 0;
		float sum_non0 = 0;

		/*
		for (unsigned int i = 0; i < n; i++) {
			cerr << cov_distr[i] << " ";
		}
		cerr << endl;
		*/
		
		sort (cov_distr.begin(), cov_distr.end());  
		
		for (unsigned int i = 0; i < n; i++) {
			sum += cov_distr[i];
			if(cov_distr[i] != 0) { sum_non0 += cov_distr[i]; n_non0++; }
			if(cov_distr[i] < min_cov) { min_cov = cov_distr[i]; }
			if(cov_distr[i] < min_nonzero_cov && cov_distr[i] != 0) { min_nonzero_cov = cov_distr[i]; }
		}
		
		if(n>0) { mean_cov =(float)sum/(float)n; }
		else { mean_cov = 0; }
		if (n_non0 > 0) { mean_cov_non0 =(float)sum_non0/(float)n_non0; }
		else { mean_cov_non0 = 0; }
		
		median_cov = cov_distr[int(ceil(n/2))];
	}
	
	void addCovNfwd(int c) { cov_distr_N_fwd.push_back(c); }	
	void addCovTfwd(int c) { cov_distr_T_fwd.push_back(c); }
	void addCovNrev(int c) { cov_distr_N_rev.push_back(c); }	
	void addCovTrev(int c) { cov_distr_T_rev.push_back(c); }
	
	void addRefCovN(int c) { ref_cov_distr_N.push_back(c); }	
	void addRefCovT(int c) { ref_cov_distr_T.push_back(c); }

	int getAvgCovNfwd() { return ceil(mean_cov_N_fwd); }
	int getAvgCovNrev() { return ceil(mean_cov_N_rev); }
	
	int getAvgCovTfwd() { return ceil(mean_cov_T_fwd); }
	int getAvgCovTrev() { return ceil(mean_cov_T_rev); }
	
	int getAvgNon0CovNfwd() { return ceil(mean_cov_N_non0_fwd); }
	int getAvgNon0CovNrev() { return ceil(mean_cov_N_non0_rev); }

	int getAvgNon0CovTfwd() { return ceil(mean_cov_T_non0_fwd); }
	int getAvgNon0CovTrev() { return ceil(mean_cov_T_non0_rev); }

	int getMinCovNfwd() { return min_cov_N_fwd; }
	int getMinCovNrev() { return min_cov_N_rev; }
	
	int getMedianCovNfwd() { return median_cov_N_fwd; }
	int getMedianCovNrev() { return median_cov_N_rev; }
	
	int getMinCovTfwd() { return min_cov_T_fwd; }
	int getMinCovTrev() { return min_cov_T_rev; }

	int getMedianCovTfwd() { return median_cov_T_fwd; }
	int getMedianCovTrev() { return median_cov_T_rev; }
	
	int getMinNon0CovNfwd() { return min_nonzero_cov_N_fwd; }
	int getMinNon0CovNrev() { return min_nonzero_cov_N_rev; }

	int getMinNon0CovTfwd() { return min_nonzero_cov_T_fwd; }
	int getMinNon0CovTrev() { return min_nonzero_cov_T_rev; }

	int getMinRefCovN() { return min_ref_cov_N; }
	int getMinRefCovT() { return min_ref_cov_T; }
	
	int getMedianRefCovN() { return median_ref_cov_N; }
	int getMedianRefCovT() { return median_ref_cov_T; }
	
	int getMinNon0RefCovN() { return min_nonzero_ref_cov_N; }
	int getMinNon0RefCovT() { return min_nonzero_ref_cov_T; }
	
	int getAvgRefCovN() { return ceil(mean_ref_cov_N); }
	int getAvgRefCovT() { return ceil(mean_ref_cov_T); }
	int getAvgNon0RefCovN() { return ceil(mean_ref_cov_N_non0); }
	int getAvgNon0RefCovT() { return ceil(mean_ref_cov_T_non0); }
};

#endif
