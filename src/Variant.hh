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

using namespace std;

class Variant_t
{
public:

	string chr;
	int pos;
	string ref;
	string alt;
	int ref_cov_normal;
	int ref_cov_tumor;
	int alt_cov_normal;
	int alt_cov_tumor;
	char prev_bp_ref; // base-pair preceding the mutation in reference
	char prev_bp_alt; // base-pair preceding the mutation in alternative


	Variant_t(string chr_, int pos_, string ref_, string alt_, int ref_cov_normal_, int ref_cov_tumor_, int alt_cov_normal_, int alt_cov_tumor_, char prev_bp_ref_, char prev_bp_alt_)
	{ 
		chr = chr_;
		pos = pos_;
		ref = ref_;
		alt = alt_;
		ref_cov_normal = ref_cov_normal_;
		ref_cov_tumor = ref_cov_tumor_;
		alt_cov_normal = alt_cov_normal_;
		alt_cov_tumor = alt_cov_tumor_;
		prev_bp_ref = prev_bp_ref_;
		prev_bp_alt = prev_bp_alt_;
	}
	
	void printVCF();
	string genotype(int R, int A);
};

#endif
