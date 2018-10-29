#include "Ref.hh"

/****************************************************************************
** Ref.cc
**
** Class for storing kmer information for the reference sequence
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


// allocate mmeory for data structures
void Ref_t::init() {
	
	// allocate memory for mertables
	mertable_nml = new 	unordered_map<string,cov_t>();
	mertable_tmr = new 	unordered_map<string,cov_t>();
	
	// allocate memory for coverage info
	normal_coverage = new vector<cov_t>(); // normal k-mer coverage across the reference
	tumor_coverage  = new vector<cov_t>(); // tumor k-mer coverage across the reference
	resetCoverage();
}

// index mers
void Ref_t::indexMers()
{
	if (!indexed_m)
	{
		
		assert(mertable_nml != NULL);
		assert(mertable_tmr != NULL);
		
		CanonicalMer_t cmer;

		//for (unsigned int i = 0; i < (seq.length() - K + 1); ++i)
		for (unsigned int i = 0; (i + K) < seq.length(); ++i)
		{
			//assert(i<seq.length()); // check for out of range index
			cmer.set(seq.substr(i, K));
			
			cov_t c;
		    mertable_nml->insert(std::pair<string,cov_t>(cmer.mer_m,c));
		    mertable_tmr->insert(std::pair<string,cov_t>(cmer.mer_m,c));
			
		}
		indexed_m = true;
	}
}

// return true if the input mer is found in the mertable
bool Ref_t::hasMer(const string & cmer)
{
	indexMers();
	return mertable_nml->count(cmer);
}

// updated coverage for input mer
void Ref_t::updateCoverage(const string & cmer, int cov, unsigned int strand, int sample) {
	indexMers();
	
	unordered_map<string,cov_t> * mertable = NULL;
		
	if(sample == TMR)      { mertable = mertable_tmr; }
	else if(sample == NML) { mertable = mertable_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	assert(mertable != NULL);
	//if (mertable == NULL) { cerr << "Error: null pointer to mer-table!" << endl; } 
	
	std::unordered_map<string,cov_t>::iterator it = mertable->find(cmer);
	if (it != mertable->end()) {
		if(strand == FWD) { 
			((*it).second).fwd = cov; 
		}
		else if(strand == REV) { 
			((*it).second).rev = cov; 
		}
	}
}

// updated haplotype coverage for input mer
void Ref_t::updateHPCoverage(const string & cmer, int hp0_cov, int hp1_cov, int hp2_cov, int sample) {
	indexMers();
	
	unordered_map<string,cov_t> * mertable = NULL;
		
	if(sample == TMR)      { mertable = mertable_tmr; }
	else if(sample == NML) { mertable = mertable_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	assert(mertable != NULL);
	//if (mertable == NULL) { cerr << "Error: null pointer to mer-table!" << endl; } 
	
	std::unordered_map<string,cov_t>::iterator it = mertable->find(cmer);
	if (it != mertable->end()) {
		((*it).second).hp0 = hp0_cov; 
		((*it).second).hp1 = hp1_cov; 
		((*it).second).hp2 = hp2_cov; 
	}
}

// compute kmer coverage over the reference sequence 
void Ref_t::computeCoverage(int sample) {
	CanonicalMer_t cmer;
	
	unordered_map<string,cov_t> * mertable = NULL;
	vector<cov_t> * coverage = NULL;
	
	if(sample == TMR)      { mertable = mertable_tmr; coverage = tumor_coverage; }
	else if(sample == NML) { mertable = mertable_nml; coverage = normal_coverage; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	assert(mertable != NULL);
	assert(coverage != NULL);

	for (unsigned i = 0; (i + K) < seq.length(); ++i) 	
	{	
		cmer.set(seq.substr(i, K));			
		std::unordered_map<string,cov_t>::iterator it = mertable->find(cmer.mer_m);
		if (it != mertable->end()) {
			int cov_fwd = ((*it).second).fwd;
			int cov_rev = ((*it).second).rev;
			int cov_hp0 = ((*it).second).hp0;
			int cov_hp1 = ((*it).second).hp1;
			int cov_hp2 = ((*it).second).hp2;
			
			/*
			normal_coverage.at(i) = n_cov;			
			tumor_coverage.at(i) = t_cov;
			
			if( i==(end-1) ) {
				for (int l = 0; l < K; l++) {
					normal_coverage.at(i+l) = n_cov;
					tumor_coverage.at(i+l) = t_cov;
				}
			}
			*/
			
			if(i==0) {
				for (int j=i; j<K; ++j) { 
					coverage->at(j).fwd = cov_fwd; 				
					coverage->at(j).fwd = cov_rev; 
					coverage->at(j).hp0 = cov_hp0;
					coverage->at(j).hp1 = cov_hp1; 				
					coverage->at(j).hp2 = cov_hp2; 
				}
			}
			else {
				coverage->at(i+K-1).fwd = cov_fwd;
				coverage->at(i+K-1).rev = cov_rev;
				coverage->at(i+K-1).hp0 = cov_hp0;
				coverage->at(i+K-1).hp1 = cov_hp1;
				coverage->at(i+K-1).hp2 = cov_hp2;
				
				//for (int l = 0; l < K-1; l++) {
					//if(normal_coverage.at(i+l) < n_cov) { normal_coverage.at(i+l) = n_cov; }
					//if(tumor_coverage.at(i+l) < t_cov) { tumor_coverage.at(i+l) = t_cov; }
					//normal_coverage.at(i+l) = n_cov;
					//tumor_coverage.at(i+l) = t_cov;
				//}
			}
			
		}
		else {
			coverage->at(i).fwd = 0;
			coverage->at(i).rev = 0;
			coverage->at(i).hp0 = 0;
			coverage->at(i).hp1 = 0;
			coverage->at(i).hp2 = 0;
		}
	}
}

// return k-mer coverage at position 
int Ref_t::getCovAt(unsigned pos, unsigned int strand, int sample) {
	
	vector<cov_t> * coverage = NULL;
	if(sample == NML) { coverage = normal_coverage; }
	else if(sample == TMR) { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }
	
	assert(coverage != NULL);
	//if (coverage == NULL) { cerr << "Error: null pointer to coverage vector!" << endl; } 
	
	int c = 0;
	if(coverage->size()>pos) {	
		if(strand == FWD) { c = coverage->at(pos).fwd; }
		else if(strand == REV) { c = coverage->at(pos).rev; }
	}
	else { c = 0; }
	
	return c;
}

// return k-mer haplotype coverage at position 
int Ref_t::getHPCovAt(unsigned pos, unsigned int hp, int sample) {
	
	vector<cov_t> * coverage = NULL;
	if(sample == NML) { coverage = normal_coverage; }
	else if(sample == TMR) { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }
	
	assert(coverage != NULL);
	//if (coverage == NULL) { cerr << "Error: null pointer to coverage vector!" << endl; } 
	
	int c = 0;
	if(coverage->size()>pos) {	
		if(hp == 0) { c = coverage->at(pos).hp0; }
		if(hp == 1) { c = coverage->at(pos).hp1; }
		if(hp == 2) { c = coverage->at(pos).hp2; }
	}
	else { c = 0; }
	
	return c;
}

// return min k-mer coverage in radius of size K bp
int Ref_t::getMinCovInKbp(unsigned pos, int K, int sample) {
		
	vector<cov_t> * cov = NULL;
	if (sample == TMR) { cov = tumor_coverage; }
	else if (sample == NML) { cov = normal_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }	
	
	assert(cov != NULL);
	//if (cov == NULL) { cerr << "Error: null pointer to coverage vector!" << endl; } 
	
	int min = 1000000000;
	if(cov->size()>=pos) {	
		for(int i=0; i<K; ++i) {
			int C = cov->at(pos+i).fwd + cov->at(pos+i).rev;
			if(C<min) { min = C; } 
		}
	}
	return min;
}

// print k-mer coverage along the reference
void Ref_t::printKmerCoverage(int sample) {

	vector<cov_t> * coverage = NULL;
	if(sample == NML) { coverage = normal_coverage; }
	else if(sample == TMR) { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; return; }
	
	assert(coverage != NULL);
	//if (coverage == NULL) { cerr << "Error: null pointer to coverage vector!" << endl; } 
	
    cerr << "cov " << sample << "+: ";
	for (unsigned i=0; i<coverage->size(); ++i) {
	    cerr << " " << coverage->at(i).fwd;
	}
	cerr << endl;
	
    cerr << "cov " << sample << "-: ";
	for (unsigned i=0; i<coverage->size(); ++i) {
	    cerr << " " << coverage->at(i).rev;
	}
	cerr << " len: " << coverage->size() << endl;
}

// clear DT and free memory
void Ref_t::clear() {
	
	if(mertable_nml != NULL)    { mertable_nml->clear();    delete mertable_nml;    mertable_nml = NULL;    }
	if(mertable_tmr != NULL)    { mertable_tmr->clear();    delete mertable_tmr;    mertable_tmr = NULL;    }
	if(normal_coverage != NULL) { normal_coverage->clear(); delete normal_coverage; normal_coverage = NULL; }
	if(tumor_coverage != NULL)  { tumor_coverage->clear();  delete tumor_coverage;  tumor_coverage = NULL;  }
}

// reset coverage to 0
void Ref_t::resetCoverage() {
	
	assert(normal_coverage != NULL);
	assert(tumor_coverage != NULL);
	//if (normal_coverage == NULL) { cerr << "Error: null pointer to coverage vector for normal!" << endl; } 
	//if (tumor_coverage == NULL) { cerr << "Error: null pointer to coverage vector for tumor!" << endl; } 
	
	normal_coverage->resize(seq.size()); 
	tumor_coverage->resize(seq.size());
	
	for (unsigned i=0; i<normal_coverage->size(); ++i) { 
		normal_coverage->at(i).fwd = 0;
		normal_coverage->at(i).rev = 0;
		tumor_coverage->at(i).fwd = 0;
		tumor_coverage->at(i).rev = 0;
	}
}
