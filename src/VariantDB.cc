#include "VariantDB.hh"

/******************************************************************
** VariantDB.cc
**
** Class for storing multiple variants (DB-style)
**
**  Authors: Giuseppe Narzisi
**    Date: December 21, 2015
**
*******************************************************************/

// add variant to DB
void VariantDB_t::addVar(Variant_t v) {
	
	string key = sha256(v.getSignature());
    map<string,Variant_t>::iterator it = DB.find(key);
	
	if (it != DB.end()) {
		Variant_t old_v = it->second;
		
		// keep highest supporting coverage found
		if (it->second.ref_cov_normal < v.ref_cov_normal) { it->second.ref_cov_normal = v.ref_cov_normal; }
		if (it->second.ref_cov_tumor  < v.ref_cov_tumor ) { it->second.ref_cov_tumor = v.ref_cov_tumor;   }
		if (it->second.alt_cov_normal < v.alt_cov_normal) { it->second.alt_cov_normal = v.alt_cov_normal; }
		if (it->second.alt_cov_tumor  < v.alt_cov_tumor ) { it->second.alt_cov_tumor = v.alt_cov_tumor;   }

		// re-gentype and score
		it->second.update();
	}
	else { 
		DB.insert(pair<string,Variant_t>(key,v));
	}
}

void VariantDB_t::printHeader(const string version, const string reference, char * date, Filters &fs) {
	
	cout << "##fileformat=VCFv4.1\n"
			"##fileDate=" << date << ""
			"##source=lancet " << version << "\n"
			"##reference=" << reference << "\n"
			"##INFO=<ID=FETS,Number=1,Type=Float,Description=\"phred-scaled p-value from the Fisher's exact test for tumor-normal allele counts\">\n"
			"##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic mutation\">\n"
			"##INFO=<ID=SHARED,Number=0,Type=Flag,Description=\"Shared mutation betweem tumor and normal\">\n"
			"##INFO=<ID=NORMAL,Number=0,Type=Flag,Description=\"Mutation present only in the normal\">\n"
			"##FILTER=<ID=MS,Description=\"Microsatellite mutation (format: #LEN#MOTIF)\">\n"
			"##FILTER=<ID=LowCovNormal,Description=\"low coverage in the normal (<" << fs.minCovNormal << ")\">\n"
			"##FILTER=<ID=HighCovNormal,Description=\"high coverage in the normal (>" << fs.maxCovNormal << ")\">\n"
			"##FILTER=<ID=LowCovTumor,Description=\"low coverage in the tumor (<" << fs.minCovTumor << ")\">\n"
			"##FILTER=<ID=HighCovTumor,Description=\"high coverage in the tumor (>" << fs.maxCovTumor << ")\">"
			"##FILTER=<ID=LowVafTumor,Description=\"low variant allele frequency in the tumor (<" << fs.minVafTumor << ")\">\n"
			"##FILTER=<ID=HighVafNormal,Description=\"high variant allele frequency in the normal (>" << fs.maxVafNormal << ")\">\n"
			"##FILTER=<ID=LowAltCntTumor,Description=\"low alternative allele count in the tumor (<" << fs.minAltCntTumor << ")\">\n"
			"##FILTER=<ID=HighAltCntNormal,Description=\"high alternative allele count in the normal (>" << fs.maxAltCntNormal << ")\">\n"
			"##FILTER=<ID=LowFisherScore,Description=\"low Fisher's exact test score for tumor-normal allele counts (<" << fs.minPhredFisher << ")\">\n"
			"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
			"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
			"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"depth supporting reference/indel at the site\">\n"
			"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Normal4 Tumor4\n";
}

// print variant in VCF format
void VariantDB_t::printToVCF(const string version, const string reference, char * date, Filters &fs) {
	
	printHeader(version,reference,date,fs);
	
	// dump map content to vector for custom sorting
	vector< pair<string,Variant_t> > myVec(DB.begin(), DB.end());
	// sort based on chromosome location
	sort(myVec.begin(),myVec.end(),byPos());

	vector< pair<string,Variant_t> >::iterator it;
	for (it=myVec.begin(); it!=myVec.end(); ++it) {
		//cerr << it->first << "\t";
		it->second.printVCF();
	}	
}
