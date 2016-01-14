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
		if (old_v.ref_cov_normal < v.ref_cov_normal) { old_v.ref_cov_normal = v.ref_cov_normal; }
		if (old_v.ref_cov_tumor < v.ref_cov_tumor) { old_v.ref_cov_tumor = v.ref_cov_tumor; }
		if (old_v.alt_cov_normal < v.alt_cov_normal) { old_v.alt_cov_normal = v.alt_cov_normal; }
		if (old_v.alt_cov_tumor < v.alt_cov_tumor) { old_v.alt_cov_tumor = v.alt_cov_tumor; }

		// re-gentype and score
		old_v.update();
		// update variant in DB
		it->second = old_v;
	}
	else { 
		DB.insert(pair<string,Variant_t>(key,v));
	}
}

// print variant in VCF format
void VariantDB_t::printToVCF() {
	
	// dump map content to vector for custom sorting
	vector< pair<string,Variant_t> > myVec(DB.begin(), DB.end());
	// sort based on chromosome location
	sort(myVec.begin(),myVec.end(),byPos());

	vector< pair<string,Variant_t> >::iterator it;
	for (it=myVec.begin(); it!=myVec.end(); ++it) {
		//cerr << it->first << "\t";
		it->second.printVCF();
	}	
	/*
	std::map<string,Variant_t>::iterator it;
	for (it=DB.begin(); it!=DB.end(); ++it) {
		//cerr << it->first << "\t";
		it->second.printVCF();
	}
	*/
}
