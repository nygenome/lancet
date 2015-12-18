#include "Variant.hh"

/******************************************************************
** Variant.cc
**
** Class for storing basic variant information storage
**
**  Authors: Giuseppe Narzisi
**    Date: November 5, 2015
**
*******************************************************************/

void Variant_t::printVCF() {
	//CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Pat4-FF-Normal-DNA      Pat4-FF-Tumor-DNA
	string ID = ".";
	float QUAL = 0.0;
	
	FET_t fet;
	int n11 = ref_cov_normal;
	int n12 = ref_cov_tumor;
	int n21 = alt_cov_normal;
	int n22 = alt_cov_tumor;
	double left = 0;
	double right = 0;
	double twotail = 0;
	
	double prob = fet.kt_fisher_exact(n11, n12, n21, n22, &left, &right, &twotail);
	double score = -10*log10(prob);
	
	string FILTER = "PASS";
	string INFO = "SOMATIC;FETS=" + dtos(score);
	string FOMRAT = "GT:AD:DP";	
	string GT_normal = genotype(ref_cov_normal,alt_cov_normal);
	string GT_tumor = genotype(ref_cov_tumor,alt_cov_tumor);
	
	string NORMAL = GT_normal + ":" + itos(ref_cov_normal) + "," + itos(alt_cov_normal) + ":" + itos(ref_cov_normal+alt_cov_normal);
	string TUMOR = GT_tumor + ":" + itos(ref_cov_tumor) + "," + itos(alt_cov_tumor) + ":" + itos(ref_cov_tumor+alt_cov_tumor);
	
	cout << chr << "\t" << pos << "\t" << ID << "\t" << ref << "\t" << alt << "\t" << QUAL << "\t" << FILTER << "\t" << INFO << "\t" << FOMRAT << "\t" << NORMAL << "\t" << TUMOR << endl;
}

// compute genotype info in VCF format (GT field)
//////////////////////////////////////////////////////////////
string Variant_t::genotype(int R, int A) {

	//0/0 - the sample is homozygous reference
	//0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
	//1/1 - the sample is homozygous alternate
	//1/2 - the sample is heterozygous, carrying 1 copy of the REF and 2 ALT alleles
	
	string GT = ""; // GT field on VCF format
	
	if(R>0 && A>0) { GT = "0/1"; }
	else if(R>0 && A==0) { GT = "0/0"; }
	else if(R==0 && A>0) { GT = "1/1"; }
	else if(R==0 && A==0) { GT = "."; }
	
	return GT;
}
