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
	DB.insert(pair<string,Variant_t>(key,v)); 
}

// print variant in VCF format
void VariantDB_t::printToVCF() {
	std::map<string,Variant_t>::iterator it;
	for (it=DB.begin(); it!=DB.end(); ++it) {
		//cerr << "key: " << it->first << endl;
		it->second.printVCF();
	}
}
