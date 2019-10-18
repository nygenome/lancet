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
#include <vector>
#include <algorithm>
#include "util.hh"

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

	cov_t min_alt_cov_N;
	cov_t min_alt_cov_T;
	
	cov_t min_non0_alt_cov_N;
	cov_t min_non0_alt_cov_T;
	
	cov_t mean_alt_cov_N;
	cov_t mean_alt_cov_T;
	
	cov_t mean_non0_alt_cov_N;
	cov_t mean_non0_alt_cov_T;
	
	cov_t min_ref_cov_N;
	cov_t min_ref_cov_T;
	
	cov_t min_non0_ref_cov_N;
	cov_t min_non0_ref_cov_T;
	
	cov_t mean_ref_cov_N;
	cov_t mean_ref_cov_T;

	cov_t mean_non0_ref_cov_N;
	cov_t mean_non0_ref_cov_T;

	vector<cov_t> alt_cov_N;
	vector<cov_t> alt_cov_T;
	vector<cov_t> ref_cov_N;
	vector<cov_t> ref_cov_T;
	
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative

	Transcript_t(int pos_, int ref_pos_, int start_pos_, char code_, char ref_, char qry_, cov_t alt_cov_nml, cov_t alt_cov_tmr, cov_t ref_cov_nml, cov_t ref_cov_tmr, char prev_bp_ref_, char prev_bp_alt_, int end_pos_, int ref_end_pos_, bool flag)
		: pos(pos_), ref_pos(ref_pos_), start_pos(start_pos_), code(code_), end_pos(end_pos_), ref_end_pos(ref_end_pos_)
	{ 	
		isSomatic = flag;
		ref = ref_;
		qry = qry_;
		alt_cov_N.push_back(alt_cov_nml);
		alt_cov_T.push_back(alt_cov_tmr);
		ref_cov_N.push_back(ref_cov_nml);
		ref_cov_T.push_back(ref_cov_tmr);
		
		min_alt_cov_N = alt_cov_nml;
		min_alt_cov_T = alt_cov_tmr;

		min_non0_alt_cov_N = min_alt_cov_N;
		min_non0_alt_cov_T = min_alt_cov_T;
				
		min_ref_cov_N = ref_cov_nml;
		min_ref_cov_T = ref_cov_tmr;
		
		min_non0_ref_cov_N = min_ref_cov_N;
		min_non0_ref_cov_T = min_ref_cov_T;
		
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	// update coverage stats for normal, tumor and the reference
	void updateStats() {				

		// normal
		computeStats(alt_cov_N, min_alt_cov_N, min_non0_alt_cov_N, mean_alt_cov_N, mean_non0_alt_cov_N);
		
		//tumor
		computeStats(alt_cov_T, min_alt_cov_T, min_non0_alt_cov_T, mean_alt_cov_T, mean_non0_alt_cov_T);
		
		// normal reference coverage
		computeStats(ref_cov_N, min_ref_cov_N, min_non0_ref_cov_N, mean_ref_cov_N, mean_non0_ref_cov_N);

		// tumor reference coverage
		computeStats(ref_cov_T, min_ref_cov_T, min_non0_ref_cov_T, mean_ref_cov_T, mean_non0_ref_cov_T);		
	}
	
	// compute coverage stats
	void computeStats(vector<cov_t> &cov_distr, cov_t &min, cov_t &min_non0, cov_t &mean, cov_t &mean_non0) {
		//float sum= 0;
		//float sum_non0 = 0;
		unsigned int n = cov_distr.size();
		//unsigned int n_non0 = 0;
		
		cov_t sum = {0,0,0,0,0,0,0};
		cov_t sum_non0 = {0,0,0,0,0,0,0};
		cov_t n_non0 = {0,0,0,0,0,0,0};

		/*
		cerr << endl << "Transcript coverage: ";		
		for (unsigned int i = 0; i < n; i++) {
			cerr << cov_distr[i] << " ";
		}
		cerr << "(node size = " << nodesize << ")" << endl;
		*/
		
		//sort (cov_distr.begin(), cov_distr.end());  
		
		for (unsigned int i = 0; i < n; ++i) {
			sum.fwd += cov_distr[i].fwd;
			sum.rev += cov_distr[i].rev;
			sum.minqv_fwd += cov_distr[i].minqv_fwd;
			sum.minqv_rev += cov_distr[i].minqv_rev;
			sum.hp0 += cov_distr[i].hp0;
			sum.hp1 += cov_distr[i].hp1;
			sum.hp2 += cov_distr[i].hp2;
			sum.hp0_minqv += cov_distr[i].hp0_minqv;
			sum.hp1_minqv += cov_distr[i].hp1_minqv;
			sum.hp2_minqv += cov_distr[i].hp2_minqv;
			
			if(cov_distr[i].fwd != 0) { sum_non0.fwd += cov_distr[i].fwd; ++n_non0.fwd; }
			if(cov_distr[i].rev != 0) { sum_non0.rev += cov_distr[i].rev; ++n_non0.rev; }
			if(cov_distr[i].minqv_fwd != 0) { sum_non0.minqv_fwd += cov_distr[i].minqv_fwd; ++n_non0.minqv_fwd; }
			if(cov_distr[i].minqv_rev != 0) { sum_non0.minqv_rev += cov_distr[i].minqv_rev; ++n_non0.minqv_rev; }
			if(cov_distr[i].hp0 != 0) { sum_non0.hp0 += cov_distr[i].hp0; ++n_non0.hp0; }
			if(cov_distr[i].hp1 != 0) { sum_non0.hp1 += cov_distr[i].hp1; ++n_non0.hp1; }
			if(cov_distr[i].hp2 != 0) { sum_non0.hp2 += cov_distr[i].hp2; ++n_non0.hp2; }
			if(cov_distr[i].hp0_minqv != 0) { sum_non0.hp0_minqv += cov_distr[i].hp0_minqv; ++n_non0.hp0_minqv; }
			if(cov_distr[i].hp1_minqv != 0) { sum_non0.hp1_minqv += cov_distr[i].hp1_minqv; ++n_non0.hp1_minqv; }
			if(cov_distr[i].hp2_minqv != 0) { sum_non0.hp2_minqv += cov_distr[i].hp2_minqv; ++n_non0.hp2_minqv; }

			if(cov_distr[i].fwd < min.fwd) { min.fwd = cov_distr[i].fwd; }
			if(cov_distr[i].rev < min.rev) { min.rev = cov_distr[i].rev; }			
			if(cov_distr[i].minqv_fwd < min.minqv_fwd) { min.minqv_fwd = cov_distr[i].minqv_fwd; }
			if(cov_distr[i].minqv_rev < min.minqv_rev) { min.minqv_rev = cov_distr[i].minqv_rev; }
			if(cov_distr[i].hp0 < min.hp0) { min.hp0 = cov_distr[i].hp0; }
			if(cov_distr[i].hp1 < min.hp1) { min.hp1 = cov_distr[i].hp1; }
			if(cov_distr[i].hp2 < min.hp2) { min.hp2 = cov_distr[i].hp2; }
			if(cov_distr[i].hp0_minqv < min.hp0_minqv) { min.hp0_minqv = cov_distr[i].hp0_minqv; }
			if(cov_distr[i].hp1_minqv < min.hp1_minqv) { min.hp1_minqv = cov_distr[i].hp1_minqv; }
			if(cov_distr[i].hp2_minqv < min.hp2_minqv) { min.hp2_minqv = cov_distr[i].hp2_minqv; }
			
			if(cov_distr[i].fwd < min_non0.fwd && cov_distr[i].fwd != 0) { min_non0.fwd = cov_distr[i].fwd; }
			if(cov_distr[i].rev < min_non0.rev && cov_distr[i].rev != 0) { min_non0.rev = cov_distr[i].rev; }
			if(cov_distr[i].minqv_fwd < min_non0.minqv_fwd && cov_distr[i].minqv_fwd != 0) { min_non0.minqv_fwd = cov_distr[i].minqv_fwd; }
			if(cov_distr[i].minqv_rev < min_non0.minqv_rev && cov_distr[i].minqv_rev != 0) { min_non0.minqv_rev = cov_distr[i].minqv_rev; }
			if(cov_distr[i].hp0 < min_non0.hp0 && cov_distr[i].hp0 != 0) { min_non0.hp0 = cov_distr[i].hp0; }
			if(cov_distr[i].hp1 < min_non0.hp1 && cov_distr[i].hp1 != 0) { min_non0.hp1 = cov_distr[i].hp1; }
			if(cov_distr[i].hp2 < min_non0.hp2 && cov_distr[i].hp2 != 0) { min_non0.hp2 = cov_distr[i].hp2; }
			if(cov_distr[i].hp0_minqv < min_non0.hp0_minqv && cov_distr[i].hp0_minqv != 0) { min_non0.hp0_minqv = cov_distr[i].hp0_minqv; }
			if(cov_distr[i].hp1_minqv < min_non0.hp1_minqv && cov_distr[i].hp1_minqv != 0) { min_non0.hp1_minqv = cov_distr[i].hp1_minqv; }
			if(cov_distr[i].hp2_minqv < min_non0.hp2_minqv && cov_distr[i].hp2_minqv != 0) { min_non0.hp2_minqv = cov_distr[i].hp2_minqv; }
		}
		
		if(n>0) { 
			mean.fwd = (float)sum.fwd/(float)n; 
			mean.rev = (float)sum.rev/(float)n; 
			mean.minqv_fwd = (float)sum.minqv_fwd/(float)n; 
			mean.minqv_rev = (float)sum.minqv_rev/(float)n; 
			mean.hp0 = (float)sum.hp0/(float)n; 
			mean.hp1 = (float)sum.hp1/(float)n; 
			mean.hp2 = (float)sum.hp2/(float)n; 
			mean.hp0_minqv = (float)sum.hp0_minqv/(float)n; 
			mean.hp1_minqv = (float)sum.hp1_minqv/(float)n; 
			mean.hp2_minqv = (float)sum.hp2_minqv/(float)n; 
		}
		else { 
			mean.fwd = 0; 		
			mean.rev = 0; 		
			mean.minqv_fwd = 0; 		
			mean.minqv_rev = 0; 		
			mean.hp0 = 0;
			mean.hp1 = 0;
			mean.hp2 = 0;
			mean.hp0_minqv = 0; 		
			mean.hp1_minqv = 0;		
			mean.hp2_minqv = 0;
		}
		
		if (n_non0.fwd > 0) { mean_non0.fwd = ceil((float)sum_non0.fwd/(float)n_non0.fwd); } else { mean_non0.fwd = 0; }
		if (n_non0.rev > 0) { mean_non0.rev = ceil((float)sum_non0.rev/(float)n_non0.rev); } else { mean_non0.rev = 0; }
		if (n_non0.minqv_fwd > 0) { mean_non0.minqv_fwd = ceil((float)sum_non0.minqv_fwd/(float)n_non0.minqv_fwd); } else { mean_non0.minqv_fwd = 0; }
		if (n_non0.minqv_rev > 0) { mean_non0.minqv_rev = ceil((float)sum_non0.minqv_rev/(float)n_non0.minqv_rev); } else { mean_non0.minqv_rev = 0; }
		if (n_non0.hp0 > 0) { mean_non0.hp0 = ceil((float)sum_non0.hp0/(float)n_non0.hp0); } else { mean_non0.hp0 = 0; }
		if (n_non0.hp1 > 0) { mean_non0.hp1 = ceil((float)sum_non0.hp1/(float)n_non0.hp1); } else { mean_non0.hp1 = 0; }
		if (n_non0.hp2 > 0) { mean_non0.hp2 = ceil((float)sum_non0.hp2/(float)n_non0.hp2); } else { mean_non0.hp2 = 0; }
		if (n_non0.hp0_minqv > 0) { mean_non0.hp0_minqv = ceil((float)sum_non0.hp0_minqv/(float)n_non0.hp0_minqv); } else { mean_non0.hp0_minqv = 0; }
		if (n_non0.hp1_minqv > 0) { mean_non0.hp1_minqv = ceil((float)sum_non0.hp1_minqv/(float)n_non0.hp1_minqv); } else { mean_non0.hp1_minqv = 0; }
		if (n_non0.hp2_minqv > 0) { mean_non0.hp2_minqv = ceil((float)sum_non0.hp2_minqv/(float)n_non0.hp2_minqv); } else { mean_non0.hp2_minqv = 0; }
		
		//median_cov = cov_distr[int(ceil(n/2))];
	}
	
	void addAltCovNml(cov_t c) { alt_cov_N.push_back(c); }
	void addAltCovTmr(cov_t c) { alt_cov_T.push_back(c); }
	void addRefCovNml(cov_t c) { ref_cov_N.push_back(c); }
	void addRefCovTmr(cov_t c) { ref_cov_T.push_back(c); }

	int getAvgCovNfwd() { if(code == 'x') { return mean_alt_cov_N.minqv_fwd; } else { return mean_alt_cov_N.fwd; } }
	int getAvgCovNrev() { if(code == 'x') { return mean_alt_cov_N.minqv_rev; } else { return mean_alt_cov_N.rev; } }
	
	int getAvgCovTfwd() { if(code == 'x') { return mean_alt_cov_T.minqv_fwd; } else { return mean_alt_cov_T.fwd; } }
	int getAvgCovTrev() { if(code == 'x') { return mean_alt_cov_T.minqv_rev; } else { return mean_alt_cov_T.rev; } }
	
	int getAvgNon0CovNfwd() { if(code == 'x') { return mean_non0_alt_cov_N.minqv_fwd; } else { return mean_non0_alt_cov_N.fwd; } }
	int getAvgNon0CovNrev() { if(code == 'x') { return mean_non0_alt_cov_N.minqv_rev; } else { return mean_non0_alt_cov_N.rev; } }

	int getAvgNon0CovTfwd() { if(code == 'x') { return mean_non0_alt_cov_T.minqv_fwd; } else { return mean_non0_alt_cov_T.fwd; } }
	int getAvgNon0CovTrev() { if(code == 'x') { return mean_non0_alt_cov_T.minqv_rev; } else { return mean_non0_alt_cov_T.rev; } }

	int getMinCovNfwd() { if(code == 'x') { return min_alt_cov_N.minqv_fwd; } else { return min_alt_cov_N.fwd; } }
	int getMinCovNrev() { if(code == 'x') { return min_alt_cov_N.minqv_rev; } else { return min_alt_cov_N.rev; } }
	
	int getMinCovTfwd() { if(code == 'x') { return min_alt_cov_T.minqv_fwd; } else { return min_alt_cov_T.fwd; } }
	int getMinCovTrev() { if(code == 'x') { return min_alt_cov_T.minqv_rev; } else { return min_alt_cov_T.rev; } }
	
	int getMinNon0CovNfwd() { if(code == 'x') { return min_non0_alt_cov_N.minqv_fwd; } else { return min_non0_alt_cov_N.fwd; } }
	int getMinNon0CovNrev() { if(code == 'x') { return min_non0_alt_cov_N.minqv_rev; } else { return min_non0_alt_cov_N.rev; } }

	int getMinNon0CovTfwd() { if(code == 'x') { return min_non0_alt_cov_T.minqv_fwd; } else { return min_non0_alt_cov_T.fwd; } }
	int getMinNon0CovTrev() { if(code == 'x') { return min_non0_alt_cov_T.minqv_rev; } else { return min_non0_alt_cov_T.rev; } }

	int getMinRefCovNfwd() { return min_ref_cov_N.fwd; }
	int getMinRefCovNrev() { return min_ref_cov_N.rev; }

	int getMinRefCovTfwd() { return min_ref_cov_T.fwd; }
	int getMinRefCovTrev() { return min_ref_cov_T.rev; }
	
	int getMinNon0RefCovNfwd() { return min_non0_ref_cov_N.fwd; }
	int getMinNon0RefCovNrev() { return min_non0_ref_cov_N.rev; }

	int getMinNon0RefCovTfwd() { return min_non0_ref_cov_T.fwd; }
	int getMinNon0RefCovTrev() { return min_non0_ref_cov_T.rev; }
	
	int getAvgRefCovNfwd() { return mean_ref_cov_N.fwd; }
	int getAvgRefCovNrev() { return mean_ref_cov_N.rev; }

	int getAvgRefCovTfwd() { return mean_ref_cov_T.fwd; }
	int getAvgRefCovTrev() { return mean_ref_cov_T.rev; }

	int getAvgNon0RefCovNfwd() { return mean_non0_ref_cov_N.fwd; }
	int getAvgNon0RefCovNrev() { return mean_non0_ref_cov_N.rev; }

	int getAvgNon0RefCovTfwd() { return mean_non0_ref_cov_T.fwd; }
	int getAvgNon0RefCovTrev() { return mean_non0_ref_cov_T.rev; }
	
	/*-----------------------------*/
	
	int getMinRefCovNhp0() { return min_ref_cov_N.hp0; }
	int getMinRefCovNhp1() { return min_ref_cov_N.hp1; }
	int getMinRefCovNhp2() { return min_ref_cov_N.hp2; }
	
	int getMinRefCovThp0() { return min_ref_cov_T.hp0; }
	int getMinRefCovThp1() { return min_ref_cov_T.hp1; }
	int getMinRefCovThp2() { return min_ref_cov_T.hp2; }
	
	int getMinCovNhp0() { if(code == 'x') {return min_alt_cov_N.hp0_minqv;} else {return min_alt_cov_N.hp0;} }
	int getMinCovNhp1() { if(code == 'x') {return min_alt_cov_N.hp1_minqv;} else {return min_alt_cov_N.hp1;} }
	int getMinCovNhp2() { if(code == 'x') {return min_alt_cov_N.hp2_minqv;} else {return min_alt_cov_N.hp2;} }
	
	int getMinCovThp0() { if(code == 'x') {return min_alt_cov_T.hp0_minqv;} else {return min_alt_cov_T.hp0;} }
	int getMinCovThp1() { if(code == 'x') {return min_alt_cov_T.hp1_minqv;} else {return min_alt_cov_T.hp1;} }
	int getMinCovThp2() { if(code == 'x') {return min_alt_cov_T.hp2_minqv;} else {return min_alt_cov_T.hp2;} }
	
	int getAvgRefCovNhp0() { return mean_ref_cov_N.hp0; }
	int getAvgRefCovNhp1() { return mean_ref_cov_N.hp1; }
	int getAvgRefCovNhp2() { return mean_ref_cov_N.hp2; }
	
	int getAvgRefCovThp0() { return mean_ref_cov_T.hp0; }
	int getAvgRefCovThp1() { return mean_ref_cov_T.hp1; }
	int getAvgRefCovThp2() { return mean_ref_cov_T.hp2; }
	
	int getAvgCovNhp0() { if(code == 'x') {return mean_alt_cov_N.hp0_minqv; } else {return mean_alt_cov_N.hp0;} }
	int getAvgCovNhp1() { if(code == 'x') {return mean_alt_cov_N.hp1_minqv; } else {return mean_alt_cov_N.hp1;} }
	int getAvgCovNhp2() { if(code == 'x') {return mean_alt_cov_N.hp2_minqv; } else {return mean_alt_cov_N.hp2;} }
	
	int getAvgCovThp0() { if(code == 'x') {return mean_alt_cov_T.hp0_minqv; } else {return mean_alt_cov_T.hp0;} }
	int getAvgCovThp1() { if(code == 'x') {return mean_alt_cov_T.hp1_minqv; } else {return mean_alt_cov_T.hp1;} }
	int getAvgCovThp2() { if(code == 'x') {return mean_alt_cov_T.hp2_minqv; } else {return mean_alt_cov_T.hp2;} }
	
};

#endif
