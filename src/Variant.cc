#include "Variant.hh"

/****************************************************************************
** Variant.cc
**
** Class for storing basic variant information storage
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


void Variant_t::printVCF() {
	//CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Pat4-FF-Normal-DNA      Pat4-FF-Tumor-DNA
	string ID = ".";
	string FILTER = "";
	
	int tot_ref_cov_tumor = ref_cov_tumor_fwd + ref_cov_tumor_rev;	
	int tot_alt_cov_tumor = alt_cov_tumor_fwd + alt_cov_tumor_rev;

	int tot_ref_cov_normal = ref_cov_normal_fwd + ref_cov_normal_rev;	
	int tot_alt_cov_normal = alt_cov_normal_fwd + alt_cov_normal_rev;
		
	double fet_score = compute_FET_score(); // compute FET score
	double fet_score_strand_bias = compute_SB_score(); // compute strand bias score
	
	double fet_score_haplotype_normal = 0.0;
	double fet_score_haplotype_tumor = 0.0;
	if (LR_MODE) {
		fet_score_haplotype_normal = compute_HP_score(HPRN[0], HPRN[1], HPAN[0], HPAN[1]);
		fet_score_haplotype_tumor = compute_HP_score(HPRT[0], HPRT[1], HPAT[0], HPAT[1]);
	}
	
	string status = "?";
	char flag = bestState(tot_ref_cov_normal,tot_alt_cov_normal,tot_ref_cov_tumor,tot_alt_cov_tumor);
	if(flag == 'T') { status = "SOMATIC"; }
	else if(flag == 'S') { status = "SHARED"; }
	else if(flag == 'L') { status = "LOH"; }
	else if(flag == 'N') { status = "NORMAL"; }
	else if(flag == 'E') { status = "NONE"; return; } // do not print varaints without support
		
	string INFO = status + ";FETS=" + dtos(fet_score);
	if(type=='I') { INFO += ";TYPE=ins"; }
	if(type=='D') { INFO += ";TYPE=del"; }
	if(type=='S') { INFO += ";TYPE=snv"; }
	if(type=='C') { INFO += ";TYPE=complex"; }
	
	INFO += ";LEN=" + itos(len) + ";KMERSIZE=" + itos(kmer) + ";SB=" + dtos(fet_score_strand_bias);
	
	if (LR_MODE) { // in linked-read mode add fet-score of haplotype bias in the tumor/normal
		INFO += ";HPSN=" + dtos(fet_score_haplotype_normal) + ";HPST=" + dtos(fet_score_haplotype_tumor) ;
	}
	
	if(!str.empty()) { INFO += ";MS=" + str; } // add STR info
	
	double QUAL = fet_score;	
	// apply filters
	
	int tumor_cov = tot_ref_cov_tumor + tot_alt_cov_tumor;
	double tumor_vaf = (tumor_cov == 0) ? 0 : ((double)tot_alt_cov_tumor/(double)tumor_cov);
		
	int normal_cov = tot_ref_cov_normal + tot_alt_cov_normal;
	double normal_vaf = (normal_cov == 0) ? 0 : ((double)tot_alt_cov_normal/(double)normal_cov);
		
	if (filters == NULL) { cerr << "Error: filters not assigned" << endl; }
	
	if(!str.empty()) { // STR variant
		if(fet_score < filters->minPhredFisherSTR) {
			if (FILTER.compare("") == 0) { FILTER = "LowFisherSTR"; }
			else { FILTER += ";LowFisherSTR"; }
		}
	}
	// non-STR variant
	else if(fet_score < filters->minPhredFisher) { 
		if (FILTER.compare("") == 0) { FILTER = "LowFisherScore"; }
		else { FILTER += ";LowFisherScore"; }
	}
	
	// the following filters are applied to all variants
	//if(fet_score < filters->minPhredFisher) { 
	//	if (FILTER.compare("") == 0) { FILTER = "LowFisherScore"; }
	//	else { FILTER += ";LowFisherScore"; }
	//}
	if(normal_cov < filters->minCovNormal) { 
		if (FILTER.compare("") == 0) { FILTER = "LowCovNormal"; }
		else { FILTER += ";LowCovNormal"; }	
	}
	if(normal_cov > filters->maxCovNormal) { 
		if (FILTER.compare("") == 0) { FILTER = "HighCovNormal"; }
		else { FILTER += ";HighCovNormal"; }	
	}
	if(tumor_cov < filters->minCovTumor) { 
		if (FILTER.compare("") == 0) { FILTER = "LowCovTumor"; }
		else { FILTER += ";LowCovTumor"; }	
	}
	if(tumor_cov > filters->maxCovTumor) { 
		if (FILTER.compare("") == 0) { FILTER = "HighCovTumor"; }
		else { FILTER += ";HighCovTumor"; }	
	}
	if(tumor_vaf < filters->minVafTumor) { 
		if (FILTER.compare("") == 0) { FILTER = "LowVafTumor"; }
		else { FILTER += ";LowVafTumor"; }	
	}
	if(normal_vaf > filters->maxVafNormal) { 
		if (FILTER.compare("") == 0) { FILTER = "HighVafNormal"; }
		else { FILTER += ";HighVafNormal"; }	
	}
	if(tot_alt_cov_tumor < filters->minAltCntTumor) { 
		if (FILTER.compare("") == 0) { FILTER = "LowAltCntTumor"; }
		else { FILTER += ";LowAltCntTumor"; }	
	}
	if(tot_alt_cov_normal > filters->maxAltCntNormal) { 
		if (FILTER.compare("") == 0) { FILTER = "HighAltCntNormal"; }
		else { FILTER += ";HighAltCntNormal"; }
	}
	
	// if only 2 reads supporting the variant check for strand bias
	/*
	if(tot_alt_cov_tumor == 2) {
		if( (alt_cov_tumor_fwd == 0) || (alt_cov_tumor_rev == 0) ) { 
			if (FILTER.compare("") == 0) { FILTER = "StrandBias"; }
			else { FILTER += ";StrandBias"; }
		}
	}
	*/
	
	// snv specific filters
	//if( (type == 'S') && (tot_alt_cov_tumor > 2) ) { // for snv strand bias filter is applied at all coverages
		if( (alt_cov_tumor_fwd < filters->minStrandBias) || (alt_cov_tumor_rev < filters->minStrandBias) ) { 
			if (FILTER.compare("") == 0) { FILTER = "StrandBias"; }
			else { FILTER += ";StrandBias"; }
		}
	//}
	
	/*
	if(!str.empty()) { 
		if (FILTER.compare("") == 0) { FILTER = "STR"; }
		else { FILTER += ";STR"; }
	}	
	*/
		
	if(FILTER.compare("") == 0) { FILTER = "PASS"; }
		
	//compute genotype	
	GT_normal = genotype((ref_cov_normal_fwd+ref_cov_normal_rev),(alt_cov_normal_fwd+alt_cov_normal_rev));
	GT_tumor = genotype((ref_cov_tumor_fwd+ref_cov_tumor_rev),(alt_cov_tumor_fwd+alt_cov_tumor_rev));
		
	string NORMAL = "";
	string TUMOR = "";
	string FORMAT = "GT:AD:SR:SA:DP";
	
	NORMAL = GT_normal + ":" + 
		itos(tot_ref_cov_normal) + "," + itos(tot_alt_cov_normal) + ":" + 
		itos(ref_cov_normal_fwd) + "," + itos(ref_cov_normal_rev) + ":" + 
		itos(alt_cov_normal_fwd) + "," + itos(alt_cov_normal_rev) + ":" + 
		itos(tot_ref_cov_normal+tot_alt_cov_normal);

	TUMOR = GT_tumor + ":" + 
		itos(tot_ref_cov_tumor) + "," + itos(tot_alt_cov_tumor) + ":" + 
		itos(ref_cov_tumor_fwd) + "," + itos(ref_cov_tumor_rev) + ":" + 
		itos(alt_cov_tumor_fwd) + "," + itos(alt_cov_tumor_rev) + ":" + 
		itos(tot_ref_cov_tumor+tot_alt_cov_tumor);

	if (LR_MODE) {
		
		FORMAT += ":HPR:HPA";
		
		string HPrefN = itos(HPRN[0]) + "," + itos(HPRN[1]) + "," + itos(HPRN[2]);
		string HPaltN = itos(HPAN[0]) + "," + itos(HPAN[1]) + "," + itos(HPAN[2]);
		string HPrefT = itos(HPRT[0]) + "," + itos(HPRT[1]) + "," + itos(HPRT[2]);
		string HPaltT = itos(HPAT[0]) + "," + itos(HPAT[1]) + "," + itos(HPAT[2]);
	
		NORMAL += ":" + HPrefN + ":" + HPaltN;
		TUMOR += ":" + HPrefT + ":" + HPaltT;
	}
	
	cout << chr << "\t" << pos << "\t" << ID << "\t" << ref << "\t" << alt << "\t" << QUAL << "\t" << FILTER << "\t" << INFO << "\t" << FORMAT << "\t" << NORMAL << "\t" << TUMOR << endl;	
}

// compute genotype info in VCF format (GT field)
//////////////////////////////////////////////////////////////
string Variant_t::genotype(int R, int A) {

	//0/0 - the sample is homozygous reference
	//0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
	//1/1 - the sample is homozygous alternate
	//1/2 - the sample is heterozygous, carrying 1 copy of the REF and 2 ALT alleles
	
	string GT = "."; // GT field in VCF format
	
	if(R>0 && A>0) { GT = "0/1"; }
	else if(R>0 && A==0) { GT = "0/0"; }
	else if(R==0 && A>0) { GT = "1/1"; }
	else if(R==0 && A==0) { GT = "."; }
	
	return GT;
}

// compute genotype info in VCF format for this variant
//////////////////////////////////////////////////////////////
void Variant_t::reGenotype() {
	
	//compute genotype
	GT_normal = genotype((ref_cov_normal_fwd+ref_cov_normal_rev),(alt_cov_normal_fwd+alt_cov_normal_rev));
	GT_tumor = genotype((ref_cov_tumor_fwd+ref_cov_tumor_rev),(alt_cov_tumor_fwd+alt_cov_tumor_rev));
}


// 	compute fisher exaxt test score for tumor/normal coverages
//////////////////////////////////////////////////////////////
double Variant_t::compute_FET_score() {
	
	double prob = 0.0;
	double left = 0.0;
	double right = 0.0;
	double twotail = 0.0;

	double fet_score = 0.0;	
	//double fet_score_right = 0.0;
	
	FET_t fet;
	
	int RCN = (ref_cov_normal_fwd+ref_cov_normal_rev);
	int RCT = (ref_cov_tumor_fwd+ref_cov_tumor_rev);
	int ACN = (alt_cov_normal_fwd+alt_cov_normal_rev);
	int ACT = (alt_cov_tumor_fwd+alt_cov_tumor_rev);
		
	//cout << "RCN:"  << RCN << " RCT:" << RCT << " ACN:" << ACN << " ACT:" << ACT << endl;
	
	// fisher exaxt test (FET) score for tumor/normal coverages
	prob = fet.kt_fisher_exact(RCN, RCT, ACN, ACT, &left, &right, &twotail);
		
	if(prob == 1.0) { fet_score = 0.0; }
	else if(prob == 0.0) { fet_score = -10.0*log10(1/std::numeric_limits<double>::max()); }
	else { fet_score = -10.0*log10(prob); }	
	//cerr << endl;
	//cerr.precision(100);
	//cerr << "FET score left: " << fet_score_left << " (p = " << left << ")" << endl; 
	//cerr << "FET score right: " << fet_score_right << " (p = " << right << ")" << endl; 
	//cerr << "FET score twotail: " << fet_score_twotail << " (p = " << twotail << ")" << endl; 
	
	return fet_score;
}

// compute fisher exact test score for haplotype bias (HP) in tumor
double Variant_t::compute_HP_score(int hpr1, int hpr2, int hpa1, int hpa2) {
	
	double prob = 0.0;
	double left = 0.0;
	double right = 0.0;
	double twotail = 0.0;
	
	double hp_score = 0.0;
	
	FET_t fet;
	
	// fisher exaxt test score for haplotype bias in tumor
	prob = fet.kt_fisher_exact(hpr1, hpr2, hpa1, hpa2, &left, &right, &twotail);
	if(prob == 1) { hp_score = 0.0; }
	else { hp_score = -10.0*log10(prob); }
	
	return hp_score;
}

// compute fisher exact test score for strand bias (SB) in tumor
double Variant_t::compute_SB_score() {
	
	double prob = 0.0;
	double left = 0.0;
	double right = 0.0;
	double twotail = 0.0;
	
	double sb_score = 0.0;
	
	FET_t fet;
	
	// fisher exaxt test score for strand bias in tumor
	prob = fet.kt_fisher_exact(ref_cov_tumor_fwd, ref_cov_tumor_rev, alt_cov_tumor_fwd, alt_cov_tumor_rev, &left, &right, &twotail);
	if(prob == 1) { sb_score = 0.0; }
	else { sb_score = -10.0*log10(prob); }
	
	return sb_score;
}

// compute best state for the variant
//////////////////////////////////////////////////////////////
char Variant_t::bestState(int Rn, int An, int Rt, int At) {
	
	char ans = '?';
	if ( (An>0) && (At>0) ) { // shared (support in both samples)
		ans = 'S';
		//if( (Rn>0) && (Rt==0) ) { ans = 'L'; } // loss of heterozygosity (LOH) [het in normal but hom in tumor]
	}
	else if ( (An==0) && (At>0) ) { // somatic (alternative support only in tumor)
		ans = 'T';
	}
	else if ( (An>0) && (At==0) ) { // normal (alternative support only in normal)
		ans = 'N';
		//if( (Rn>0) && (Rt>0) ) { ans = 'L'; } // loss of heterozygosity (LOH) [het in normal but lost in tumor]
	}
	else if ( (An==0) && (At==0) ) { // no support
		ans = 'E';
	}
	return ans;
}

string Variant_t::getSignature() const {
		
	string ans = chr+":"+itos(pos)+":"+type+":"+itos(len)+":"+ref+":"+alt;
	//cerr << ans << endl;
	return ans;
}

string Variant_t::getPosition() {
		
	string ans = chr+":"+itos(pos);
	//cerr << ans << endl;
	return ans;
}