#ifndef VARIANT_HH
#define VARIANT_HH 1

/****************************************************************************
** Variant.hh
**
** Class for storing basic variant information storage
**
**  Authors: Giuseppe Narzisi
**    Date: November 5, 2015
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
#include <iostream>
#include "util.hh"
#include "FET.hh"

using namespace std;

struct Filters
{
	// filter thresholds
	double minPhredFisher;
	double maxVafNormal;
	double minVafTumor;
	int minCovNormal;
	int maxCovNormal;
	int minCovTumor;
	int maxCovTumor;
	int minAltCntTumor;
	int maxAltCntNormal;
};

class Variant_t
{
public:

	string chr;
	int pos;
	char type;
	int len;
	string ref;
	string alt;
	char status; // T=somatic, S=shared
	int ref_cov_normal;
	int ref_cov_tumor;
	int alt_cov_normal;
	int alt_cov_tumor;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative
	string GT_normal;
	string GT_tumor;
	double fet_score;
	
	// filter thresholds
	double minPhredFisher;
	double maxVafNormal;
	int minCovNormal;
	int maxCovNormal;
	int minCovTumor;
	int maxCovTumor;
	int minVafTumor;
	int minAltCntTumor;
	int maxAltCntNormal;
	
	Filters filters; // filter thresholds

	Variant_t(string chr_, int pos_, string ref_, string alt_, int ref_cov_normal_, int ref_cov_tumor_, int alt_cov_normal_, int alt_cov_tumor_, char prev_bp_ref_, char prev_bp_alt_, Filters &fs)
	{ 	
		filters = fs;
		
		chr = chr_;
		pos = pos_;
		if(ref_.at(0) == '-') { type = 'I'; ref_ = ""; len = alt_.length(); }  // deletion
		if(alt_.at(0) == '-') { type = 'D'; alt_ = ""; len = ref_.length(); }  // insertion
		if(ref_.size()==1 && alt_.size()==1 && ref_.at(0)!='-' && alt_.at(0)!='-') { type = 'S'; pos++; } // snp 
		if(type != 'S') {
			ref = prev_bp_alt_ + ref_;
			alt = prev_bp_alt_ + alt_;
		}
		else { alt = alt_; ref = ref_; len = 1; }
		
		ref_cov_normal = ref_cov_normal_;
		ref_cov_tumor = ref_cov_tumor_;
		alt_cov_normal = alt_cov_normal_;
		alt_cov_tumor = alt_cov_tumor_;
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
		
		// update score and genotypes
		update();
	}
	
	void printVCF();
	string genotype(int R, int A);
	char bestState(int Rn, int An, int Rt, int At);
	string getSignature();
	void update();
};

#endif
