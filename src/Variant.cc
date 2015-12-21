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
	
	string FILTER = "PASS";
	string INFO = "SOMATIC;FETS=" + dtos(fet_score);
	if(type=='I') { INFO += ";TYPE=ins"; }
	if(type=='D') { INFO += ";TYPE=del"; }
	if(type=='S') { INFO += ";TYPE=snp"; }
	
	double QUAL = fet_score;
	string FORMAT = "GT:AD:DP";	
	string GT_normal = genotype(ref_cov_normal,alt_cov_normal);
	string GT_tumor = genotype(ref_cov_tumor,alt_cov_tumor);
	
	string NORMAL = GT_normal + ":" + itos(ref_cov_normal) + "," + itos(alt_cov_normal) + ":" + itos(ref_cov_normal+alt_cov_normal);
	string TUMOR = GT_tumor + ":" + itos(ref_cov_tumor) + "," + itos(alt_cov_tumor) + ":" + itos(ref_cov_tumor+alt_cov_tumor);
	
	cout << chr << "\t" << pos << "\t" << ID << "\t" << ref << "\t" << alt << "\t" << QUAL << "\t" << FILTER << "\t" << INFO << "\t" << FORMAT << "\t" << NORMAL << "\t" << TUMOR << endl;
}

// compute genotype info in VCF format (GT field)
//////////////////////////////////////////////////////////////
string Variant_t::genotype(int R, int A) {

	//0/0 - the sample is homozygous reference
	//0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
	//1/1 - the sample is homozygous alternate
	//1/2 - the sample is heterozygous, carrying 1 copy of the REF and 2 ALT alleles
	
	string GT = ""; // GT field in VCF format
	
	if(R>0 && A>0) { GT = "0/1"; }
	else if(R>0 && A==0) { GT = "0/0"; }
	else if(R==0 && A>0) { GT = "1/1"; }
	else if(R==0 && A==0) { GT = "."; }
	
	return GT;
}

string Variant_t::getSignature() {
		
	string ans = chr+":"+itos(pos)+":"+type+":"+itos(len)+":"+ref+":"+alt; 
	return ans;
}
