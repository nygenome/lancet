#ifndef VARIANT_HH
#define VARIANT_HH 1

/******************************************************************
** Variant.hh
**
** Class for storing basic variant information storage
**
**  Authors: Giuseppe Narzisi
**    Date: November 5, 2015
**
*******************************************************************/

#include <string>
#include <iostream>
#include "util.hh"
#include "FET.hh"

using namespace std;

class Variant_t
{
public:

	string chr;
	int pos;
	char type;
	int len;
	string ref;
	string alt;
	int ref_cov_normal;
	int ref_cov_tumor;
	int alt_cov_normal;
	int alt_cov_tumor;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative
	string GT_normal;
	string GT_tumor;
	double fet_score;

	Variant_t(string chr_, int pos_, string ref_, string alt_, int ref_cov_normal_, int ref_cov_tumor_, int alt_cov_normal_, int alt_cov_tumor_, char prev_bp_ref_, char prev_bp_alt_)
	{ 	
		chr = chr_;
		pos = pos_;
		if(ref_.at(0) == '-') { type = 'D'; ref_ = ""; len = alt_.length(); }  // deletion
		if(alt_.at(0) == '-') { type = 'I'; alt_ = ""; len = ref_.length(); }  // insertion
		if(ref_.size()==1 && alt_.size()==1 && ref_.at(0)!='-' && alt_.at(0)!='-') { type = 'S'; } // snp 
		if(type != 'S') {
			ref = prev_bp_alt_ + ref_;
			alt = prev_bp_alt_ + alt_;
			len = 1;
		}
		else { alt = alt_; ref = ref_; }
		
		ref_cov_normal = ref_cov_normal_;
		ref_cov_tumor = ref_cov_tumor_;
		alt_cov_normal = alt_cov_normal_;
		alt_cov_tumor = alt_cov_tumor_;
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
		
		double left = 0;
		double right = 0;
		double twotail = 0;
		FET_t fet;
		double prob = fet.kt_fisher_exact(ref_cov_normal, ref_cov_tumor, alt_cov_normal, alt_cov_tumor, &left, &right, &twotail);
		fet_score = -10*log10(prob);
		
		//compute genotype
		GT_normal = genotype(ref_cov_normal,alt_cov_normal);
		GT_tumor = genotype(ref_cov_tumor,alt_cov_tumor);
		//fet_score = Fisher_score(ref_cov_normal,ref_cov_tumor,alt_cov_normal,alt_cov_tumor); // Fisher's exact test score
	}
	
	void printVCF();
	string genotype(int R, int A);
	string getSignature();
};

#endif
