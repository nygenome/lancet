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
void VariantDB_t::addVar(Variant_t v) {
	
	string key = sha256(v.getSignature());	
    map<string,Variant_t>::iterator it_v = DB.find(key);
	
	if (it_v != DB.end()) {		
		// keep highest supporting coverage found
		if (it_v->second.ref_cov_normal_fwd < v.ref_cov_normal_fwd) { it_v->second.ref_cov_normal_fwd = v.ref_cov_normal_fwd; }
		if (it_v->second.ref_cov_normal_rev < v.ref_cov_normal_rev) { it_v->second.ref_cov_normal_rev = v.ref_cov_normal_rev; }
		if (it_v->second.ref_cov_tumor_fwd  < v.ref_cov_tumor_fwd ) { it_v->second.ref_cov_tumor_fwd  = v.ref_cov_tumor_fwd;  }
		if (it_v->second.ref_cov_tumor_rev  < v.ref_cov_tumor_rev ) { it_v->second.ref_cov_tumor_rev  = v.ref_cov_tumor_rev;  }
		if (it_v->second.alt_cov_normal_fwd < v.alt_cov_normal_fwd) { it_v->second.alt_cov_normal_fwd = v.alt_cov_normal_fwd; }
		if (it_v->second.alt_cov_normal_rev < v.alt_cov_normal_rev) { it_v->second.alt_cov_normal_rev = v.alt_cov_normal_rev; }		
		if (it_v->second.alt_cov_tumor_fwd  < v.alt_cov_tumor_fwd ) { it_v->second.alt_cov_tumor_fwd  = v.alt_cov_tumor_fwd;  }
		if (it_v->second.alt_cov_tumor_rev  < v.alt_cov_tumor_rev ) { it_v->second.alt_cov_tumor_rev  = v.alt_cov_tumor_rev;  }
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
	
	cout << "##fileformat=VCFv4.1\n"
			"##fileDate=" << date << ""
			"##source=lancet " << version << "\n"
			"##reference=" << reference << "\n"
			"##INFO=<ID=FETS,Number=1,Type=Float,Description=\"phred-scaled p-value of the Fisher's exact test for tumor-normal allele counts (right-sided)\">\n"
			"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n"
			"##INFO=<ID=SHARED,Number=0,Type=Flag,Description=\"Shared mutation betweem tumor and normal\">\n"
			"##INFO=<ID=NORMAL,Number=0,Type=Flag,Description=\"Mutation present only in the normal\">\n"
			"##INFO=<ID=NONE,Number=0,Type=Flag,Description=\"Mutation not supported by data\">\n"
			"##INFO=<ID=KMERSIZE,Number=1,Type=Integer,Description=\"K-mer size used to assemble the locus\">\n"
			"##INFO=<ID=SB,Number=1,Type=Integer,Description=\"Phred-scaled strand bias of the Fisher's exact test (two-sided)\">\n"
			"##INFO=<ID=MS,Number=1,Type=String,Description=\"Microsatellite mutation (format: #LEN#MOTIF)\">\n"
			"##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"Variant size in base pairs\">\n"
			"##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type (snv, del, ins)\">\n"
			"##FILTER=<ID=LowCovNormal,Description=\"low coverage in the normal (<" << fs.minCovNormal << ")\">\n"
			"##FILTER=<ID=HighCovNormal,Description=\"high coverage in the normal (>" << fs.maxCovNormal << ")\">\n"
			"##FILTER=<ID=LowCovTumor,Description=\"low coverage in the tumor (<" << fs.minCovTumor << ")\">\n"
			"##FILTER=<ID=HighCovTumor,Description=\"high coverage in the tumor (>" << fs.maxCovTumor << ")\">\n"
			"##FILTER=<ID=LowVafTumor,Description=\"low variant allele frequency in the tumor (<" << fs.minVafTumor << ")\">\n"
			"##FILTER=<ID=HighVafNormal,Description=\"high variant allele frequency in the normal (>" << fs.maxVafNormal << ")\">\n"
			"##FILTER=<ID=LowAltCntTumor,Description=\"low alternative allele count in the tumor (<" << fs.minAltCntTumor << ")\">\n"
			"##FILTER=<ID=HighAltCntNormal,Description=\"high alternative allele count in the normal (>" << fs.maxAltCntNormal << ")\">\n"
			"##FILTER=<ID=LowFisherScore,Description=\"low Fisher's exact test score for tumor-normal allele counts (<" << fs.minPhredFisher << ")\">\n"
			"##FILTER=<ID=LowFisherSTR,Description=\"low Fisher's exact test score for tumor-normal STR allele counts (<" << fs.minPhredFisherSTR << ")\">\n"
			"##FILTER=<ID=StrandBias,Description=\"strand bias: # of non-reference reads in either forward or reverse strand below threshold (<" << fs.minStrandBias << ")\">\n"
			"##FILTER=<ID=STR,Description=\"Microsatellite mutation\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"allele depth: # of supporting reference,mutation reads at the site\">\n"
			"##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"strand counts for ref: # of supporting forward,reverse reads for reference allele\">\n"
			"##FORMAT=<ID=SA,Number=.,Type=Integer,Description=\"strand counts for alt: # of supporting forward,reverse reads for alterantive allele\">\n"
			"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name_N << "\t" << sample_name_T << "\n";
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
