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

void Ref_t::indexMers()
{
	if (!indexed_m)
	{
		CanonicalMer_t cmer;
		
		// allocate memory for mertables
		mertable_nml = new 	map<string,cov_t>();
		mertable_tmr = new 	map<string,cov_t>();
		
		// allocate memory for coverage info
		normal_coverage = new vector<cov_t>(); // normal k-mer coverage across the reference
		tumor_coverage  = new vector<cov_t>(); // tumor k-mer coverage across the reference
		resetCoverage();
		
		for (unsigned int i = 0; i < seq.length() - K + 1; i++)
		{
			cmer.set(seq.substr(i, K));
		    //mertable.insert(std::pair<string,int>(cmer.mer_m,0));
		    //mertable.insert(std::pair<string,std::pair<int,int>>(cmer.mer_m,std::make_pair(0,0)));
			
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
void Ref_t::updateCoverage(const string & cmer, unsigned int strand, char sample) {
	indexMers();
	
	map<string,cov_t> * mertable = NULL;
		
	if(sample == 'T')      { mertable = mertable_tmr; }
	else if(sample == 'N') { mertable = mertable_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	std::map<string,cov_t>::iterator it = mertable->find(cmer);
	if (it != mertable->end()) {
		if(strand == FWD) { ((*it).second).fwd += 1; }
		if(strand == REV) { ((*it).second).rev += 1; }
		
		//if(sample == 'N') { ((*it).second).first += 1; }
		//else if(sample == 'T') { ((*it).second).second += 1; }
		//else { cerr << "Error: unknown sample " << sample << endl; return; }
	}
}

// compute kmer coverage over the reference sequence 
void Ref_t::computeCoverage(char sample) {
	CanonicalMer_t cmer;
	
	map<string,cov_t> * mertable = NULL;
	vector<cov_t> * coverage = NULL;
	
	if(sample == 'T')      { mertable = mertable_tmr; coverage = tumor_coverage; }
	else if(sample == 'N') { mertable = mertable_nml; coverage = normal_coverage; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }

	unsigned int end = seq.length() - K + 1;
	for (unsigned i = 0; i < end; i++) 
	{	
		cmer.set(seq.substr(i, K));			
		std::map<string,cov_t>::iterator it = mertable->find(cmer.mer_m);
		if (it != mertable->end()) {
			int cov_fwd = ((*it).second).fwd;
			int cov_rev = ((*it).second).rev;
			
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
				for (int j=i; j<K; j++) { 
					coverage->at(j).fwd = cov_fwd; 				
					coverage->at(j).fwd = cov_rev; 
				}
			}
			else {
				coverage->at(i+K-1).fwd = cov_fwd;
				coverage->at(i+K-1).rev = cov_rev;
				
				
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
		}
	}
}

// return k-mer coverage at position 
int Ref_t::getCovAt(unsigned pos, unsigned int strand, char sample) {
	
	vector<cov_t> * coverage = NULL;
	if(sample == 'N') { coverage = normal_coverage; }
	else if(sample == 'T') { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }
	
	int c = 0;
	if(coverage->size()>pos) {	
		if(strand == FWD) { c = coverage->at(pos).fwd; }
		if(strand == REV) { c = coverage->at(pos).rev; }
	}
	else { c = -1; }
	
	return c;
}

// return min k-mer coverage in radius of size K bp
int Ref_t::getMinCovInKbp(unsigned pos, int K, char sample) {
		
	vector<cov_t> * cov = NULL;
	if (sample == 'T') { cov = tumor_coverage; }
	else if (sample == 'N') { cov = normal_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }	
	
	int min = 1000000000;
	if(cov->size()>=pos) {	
		for(int i=0; i<K; i++) {
			int C = cov->at(pos+i).fwd + cov->at(pos+i).rev;
			if(C<min) { min = C; } 
		}
	}
	return min;
}

// print k-mer coverage along the reference
void Ref_t::printKmerCoverage(char sample) {

	vector<cov_t> * coverage = NULL;
	if(sample == 'N') { coverage = normal_coverage; }
	else if(sample == 'T') { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; return; }
	
    cerr << "cov " << sample << "+: ";
	for (unsigned i=0; i<coverage->size(); i++) {
	    cerr << " " << coverage->at(i).fwd;
	}
    cerr << "cov " << sample << "-: ";
	for (unsigned i=0; i<coverage->size(); i++) {
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
	
	normal_coverage->resize(seq.size()); 
	tumor_coverage->resize(seq.size());
	
	for (unsigned i=0; i<normal_coverage->size(); i++) { 
		normal_coverage->at(i).fwd = 0;
		normal_coverage->at(i).rev = 0;
		tumor_coverage->at(i).fwd = 0;
		tumor_coverage->at(i).rev = 0;
	}
}
