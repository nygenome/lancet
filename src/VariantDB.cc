#include <vector>

#include "VariantDB.hh"

/****************************************************************************
** VariantDB.cc
**
** Class for storing multiple variants (DB-style)
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

// add variant to DB and update counts per position
void VariantDB_t::addVar(const Variant_t & v) {
	
	string key = sha256(v.getSignature());	
    map<string,Variant_t>::iterator it_v = DB.find(key);	
	
	//bool flag = false;
	
	if (it_v != DB.end()) {
		// keep variant at location with highest total coverage (tumor + normal)
		
		int old_ref_cov_normal = it_v->second.ref_cov_normal_fwd + it_v->second.ref_cov_normal_rev;
		int old_ref_cov_tumor = it_v->second.ref_cov_tumor_fwd + it_v->second.ref_cov_tumor_rev;
		int old_alt_cov_normal = it_v->second.alt_cov_normal_fwd + it_v->second.alt_cov_normal_rev;
		int old_alt_cov_tumor = it_v->second.alt_cov_tumor_fwd + it_v->second.alt_cov_tumor_rev;
		
		int new_ref_cov_normal = v.ref_cov_normal_fwd + v.ref_cov_normal_rev;
		int new_ref_cov_tumor = v.ref_cov_tumor_fwd + v.ref_cov_tumor_rev;
		int new_alt_cov_normal = v.alt_cov_normal_fwd + v.alt_cov_normal_rev;
		int new_alt_cov_tumor = v.alt_cov_tumor_fwd + v.alt_cov_tumor_rev;
		
		int old_tot_cov = old_ref_cov_normal + old_ref_cov_tumor + old_alt_cov_normal + old_alt_cov_tumor;
		int new_tot_cov = new_ref_cov_normal + new_ref_cov_tumor + new_alt_cov_normal + new_alt_cov_tumor;
		
		if(old_tot_cov < new_tot_cov) {
			it_v->second.kmer = v.kmer;
			
			it_v->second.ref_cov_normal_fwd = v.ref_cov_normal_fwd;
			it_v->second.ref_cov_normal_rev = v.ref_cov_normal_rev;
			it_v->second.ref_cov_tumor_fwd  = v.ref_cov_tumor_fwd;
			it_v->second.ref_cov_tumor_rev  = v.ref_cov_tumor_rev;
			it_v->second.alt_cov_normal_fwd = v.alt_cov_normal_fwd; 
			it_v->second.alt_cov_normal_rev = v.alt_cov_normal_rev;
			it_v->second.alt_cov_tumor_fwd  = v.alt_cov_tumor_fwd;
			it_v->second.alt_cov_tumor_rev  = v.alt_cov_tumor_rev;
			
			it_v->second.HPRN = v.HPRN;
			it_v->second.HPRT = v.HPRT;
			it_v->second.HPAN = v.HPAN;
			it_v->second.HPAT = v.HPAT;
			
			//it_v->second.reGenotype(); // recompute genotype
		}
	}
	else { 
		DB.insert(pair<string,Variant_t>(key,v));
	}

	// updated counts of variants per position in the normal
	/*
    if( ((v.getGenotypeNormal()).compare("0/0")!=0) && ((v.ref_cov_normal_fwd + v.ref_cov_normal_rev)>0) ) {
		
		string pos = v.getPosition();
	    unordered_map<string,int>::iterator it_p = nCNT.find(pos);	
		
		if (it_p != nCNT.end()) { 
			it_p->second += 1; 
		}
		else { 
			nCNT.insert(pair<string,int>(pos,1)); 
		}
    }
	*/
}

void VariantDB_t::printHeader(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T) {
	
    stringstream hdr;
	
	hdr << "##fileformat=VCFv4.2\n"
			"##fileDate=" << date << ""
			"##source=lancet " << version << "\n"
			"##cmdline=" << command_line << "\n"
			"##reference=" << reference << "\n"
			"##INFO=<ID=FETS,Number=1,Type=Float,Description=\"Phred-scaled p-value of the Fisher's exact test for tumor-normal allele counts\">\n"
			"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n"
			"##INFO=<ID=SHARED,Number=0,Type=Flag,Description=\"Shared mutation betweem tumor and normal\">\n"
			"##INFO=<ID=NORMAL,Number=0,Type=Flag,Description=\"Mutation present only in the normal\">\n"
			"##INFO=<ID=NONE,Number=0,Type=Flag,Description=\"Mutation not supported by data\">\n"
			"##INFO=<ID=KMERSIZE,Number=1,Type=Integer,Description=\"K-mer size used to assemble the locus\">\n"
			"##INFO=<ID=SB,Number=1,Type=Float,Description=\"Strand bias score: phred-scaled p-value of the Fisher's exact test for the forward/reverse read counts in the tumor\">\n"
			"##INFO=<ID=MS,Number=1,Type=String,Description=\"Microsatellite mutation (format: #LEN#MOTIF)\">\n"
			"##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Variant size in base pairs\">\n"
			"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type (snv, del, ins, complex)\">\n";
	
	if(LR_MODE)	{
		hdr << "##FORMAT=<ID=HPSN,Number=1,Type=Float,Description=\"Normal haplotype score: phred-scaled p-value of the Fisher's exact test for ref/alt haplotype counts in the normal\">\n"
			   "##FORMAT=<ID=HPST,Number=1,Type=Float,Description=\"Tumor haplotype score: phred-scaled p-value of the Fisher's exact test for ref/alt haplotype counts in the tumor\">\n";
	}
				
	hdr <<	"##FILTER=<ID=LowCovNormal,Description=\"Low coverage in the normal (<" << fs.minCovNormal << ")\">\n"
			"##FILTER=<ID=HighCovNormal,Description=\"High coverage in the normal (>" << fs.maxCovNormal << ")\">\n"
			"##FILTER=<ID=LowCovTumor,Description=\"Low coverage in the tumor (<" << fs.minCovTumor << ")\">\n"
			"##FILTER=<ID=HighCovTumor,Description=\"High coverage in the tumor (>" << fs.maxCovTumor << ")\">\n"
			"##FILTER=<ID=LowVafTumor,Description=\"Low variant allele frequency in the tumor (<" << fs.minVafTumor << ")\">\n"
			"##FILTER=<ID=HighVafNormal,Description=\"High variant allele frequency in the normal (>" << fs.maxVafNormal << ")\">\n"
			"##FILTER=<ID=LowAltCntTumor,Description=\"Low alternative allele count in the tumor (<" << fs.minAltCntTumor << ")\">\n"
			"##FILTER=<ID=HighAltCntNormal,Description=\"High alternative allele count in the normal (>" << fs.maxAltCntNormal << ")\">\n"
			"##FILTER=<ID=LowFisherScore,Description=\"Low Fisher's exact test score for tumor-normal allele counts (<" << fs.minPhredFisher << ")\">\n"
			"##FILTER=<ID=LowFisherSTR,Description=\"Low Fisher's exact test score for tumor-normal STR allele counts (<" << fs.minPhredFisherSTR << ")\">\n"
			"##FILTER=<ID=StrandBias,Description=\"Strand bias: # of non-reference reads in either forward or reverse strand below threshold (<" << fs.minStrandBias << ")\">\n"
			"##FILTER=<ID=STR,Description=\"Microsatellite mutation\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allele depth: # of supporting ref,alt reads at the site\">\n"
			"##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Strand counts for ref: # of supporting forward,reverse reads for reference allele\">\n"
			"##FORMAT=<ID=SA,Number=.,Type=Integer,Description=\"Strand counts for alt: # of supporting forward,reverse reads for alterantive allele\">\n";

	if(LR_MODE)	{
		hdr << "##FORMAT=<ID=HPR,Number=.,Type=Integer,Description=\"Haplotype counts for ref: # of reads supporting reference allele in haplotype 1, 2, and 0 respectively (0 = unassigned)\">\n"
			   "##FORMAT=<ID=HPA,Number=.,Type=Integer,Description=\"Haplotype counts for alt: # of reads supporting alternative allele in haplotype 1, 2, and 0 respectively (0 = unassigned)\">\n";
	}
	
	hdr << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name_N << "\t" << sample_name_T << "\n";
	
	cout << hdr.str();
}

// print variant in VCF format
void VariantDB_t::printToVCF(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T) {
	
	cerr << "Export variants to VCF file" << endl;
	
	printHeader(version,reference,date,fs,sample_name_N,sample_name_T);
	
	// dump map content to vector for custom sorting
	vector< pair<string,Variant_t> > myVec(DB.begin(), DB.end());
	// sort based on chromosome location
	sort(myVec.begin(),myVec.end(),byPos());

	vector< pair<string,Variant_t> >::iterator it;
	for (it=myVec.begin(); it!=myVec.end(); ++it) {
		//cerr << it->first << "\t";
		
		//string pos = (it->second).getPosition();
	    //unordered_map<string,int>::iterator itp = nCNT.find(pos);
		//if (itp == nCNT.end()) { // print variant if no muations in the normal at locus
			it->second.printVCF();
		//}
	}	
}
