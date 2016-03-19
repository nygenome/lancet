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
		for (unsigned int i = 0; i < seq.length() - K + 1; i++)
		{
			cmer.set(seq.substr(i, K));
		    //mertable.insert(std::pair<string,int>(cmer.mer_m,0));
		    mertable.insert(std::pair<string,std::pair<int,int>>(cmer.mer_m,std::make_pair(0,0)));
		}
		indexed_m = true;
	}
}

// return true if the input mer is found in the mertable
bool Ref_t::hasMer(const string & cmer)
{
	indexMers();
	return mertable.count(cmer);
}

// updated coverage for input mer
void Ref_t::updateCoverage(const string & cmer, char sample) {
	indexMers();
	std::map<string,std::pair<int,int>>::iterator it = mertable.find(cmer);
	if (it != mertable.end()) {
		if(sample == 'N') { ((*it).second).first += 1; }
		else if(sample == 'T') { ((*it).second).second += 1; }
		else { cerr << "Error: unknown sample " << sample << endl; return; }
	}
}

// compute kmer coverage over the reference sequence 
void Ref_t::computeCoverage() {
	CanonicalMer_t cmer;

	unsigned int end = seq.length() - K + 1;
	for (unsigned i = 0; i < end; i++) 
	{	
		cmer.set(seq.substr(i, K));			
		std::map<string,std::pair<int,int>>::iterator it = mertable.find(cmer.mer_m);
		if (it != mertable.end()) {
			int n_cov = ((*it).second).first;
			int t_cov = ((*it).second).second;
			
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
					normal_coverage.at(j) = n_cov; 				
					tumor_coverage.at(j) = t_cov; 
				}
			}
			else {
				normal_coverage.at(i+K-1) = n_cov;
				tumor_coverage.at(i+K-1) = t_cov;
				
				for (int l = 0; l < K-1; l++) {
					//if(normal_coverage.at(i+l) < n_cov) { normal_coverage.at(i+l) = n_cov; }
					//if(tumor_coverage.at(i+l) < t_cov) { tumor_coverage.at(i+l) = t_cov; }
					normal_coverage.at(i+l) = n_cov;
					tumor_coverage.at(i+l) = t_cov;
				}	
			}
		}
		else {
			normal_coverage.at(i) = 0;
			tumor_coverage.at(i) = 0;
		}
	}
}

// return k-mer coverage at position 
int Ref_t::getCovAt(unsigned pos, char sample) {
	 
	int c = 0;
	if(normal_coverage.size()>=pos) {	
		if (sample == 'T') { c = tumor_coverage.at(pos); }
		else if (sample == 'N') { c = normal_coverage.at(pos); }
		else { cerr << "Error: unknown sample " << sample << endl; }
	}
	else { c = -1; }
	return c;
}

// return min k-mer coverage in radius of size K bp
int Ref_t::getMinCovInKbp(unsigned pos, int K, char sample) {
		
	vector<int> cov;
	if (sample == 'T') { cov = tumor_coverage; }
	else if (sample == 'N') { cov = normal_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; }	
	
	int min = 1000000000;
	if(cov.size()>=pos) {	
		for(int i=0; i<K; i++) {
			if(cov[pos+i]<min) { min = cov[pos+i]; } 
		}
	}
	return min;
}

// print k-mer coverage along the reference
void Ref_t::printKmerCoverage(char sample) {
	
	vector<int> coverage;
	if(sample == 'N') { coverage = normal_coverage; }
	else if(sample == 'T') { coverage = tumor_coverage; }
	else { cerr << "Error: unknown sample " << sample << endl; return; }
	
    cerr << "cov " << sample << ": ";
	for (unsigned i=0; i<coverage.size(); i++) {
	    cerr << " " << coverage.at(i);
	}
	cerr << " len: " << coverage.size() << endl;
}

// reset coverage to 0
void Ref_t::resetCoverage() {
	for (unsigned i=0; i<normal_coverage.size(); i++) { 
		normal_coverage.at(i) = 0;
		tumor_coverage.at(i) = 0;
	}
}
