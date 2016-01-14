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

	string status = "?";
	char flag = bestState(ref_cov_normal,alt_cov_normal,ref_cov_tumor,alt_cov_tumor);
	if(flag == 'T') { status = "SOMATIC"; }
	else if(flag == 'S') { status = "SHARED"; }
	else if(flag == 'N') { status = "NORMAL"; }
	
	string INFO = status + ";FETS=" + dtos(fet_score);
	if(type=='I') { INFO += ";TYPE=ins"; }
	if(type=='D') { INFO += ";TYPE=del"; }
	if(type=='S') { INFO += ";TYPE=snv"; }
	
	double QUAL = fet_score;
	string FORMAT = "GT:AD:DP";		
	if(fet_score < minPhredFisher) { FILTER = "LowFisherScore"; }
	
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

void Variant_t::update() {
	
	double left = 0;
	double right = 0;
	double twotail = 0;
	FET_t fet;
	double prob = fet.kt_fisher_exact(ref_cov_normal, ref_cov_tumor, alt_cov_normal, alt_cov_tumor, &left, &right, &twotail);
	if(prob == 1) { fet_score = 0; }
	else { fet_score = -10*log10(prob); }
	
	//compute genotype
	GT_normal = genotype(ref_cov_normal,alt_cov_normal);
	GT_tumor = genotype(ref_cov_tumor,alt_cov_tumor);
}

// compute best state for the variant
//////////////////////////////////////////////////////////////
char Variant_t::bestState(int Rn, int An, int Rt, int At) {
	
	char ans = '?';
	if ( (An>0) && (At>0) ) { // shred
		ans = 'S';
	}
	else if ( (An==0) && (At>0) ) { // somatic
		ans = 'T';
	}
	else if ( (An>0) && (At==0) ) { // somatic
		ans = 'N';
	}
	return ans;
}

string Variant_t::getSignature() {
		
	string ans = chr+":"+itos(pos)+":"+type+":"+itos(len)+":"+ref+":"+alt;
	//cerr << ans << endl;
	return ans;
}
