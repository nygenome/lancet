#ifndef ERRORCORRECTOR_HH
#define ERRORCORRECTOR_HH 1

/****************************************************************************
** Variant.hh
**
** Class for error correcting low quality bases in the reads
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

//#include <string>
#include <iostream>
#include <unordered_map>

using namespace std;

class ErrorCorrector
{
public:

	ErrorCorrector() { }
	
	void mersRecovery(MerTable_t & nodes_m, int MIN_SUPPORT, int MIN_QV) {
		
		//cerr << "mers recovery" << endl;
	
		char BP[4] = { 'A', 'C', 'G', 'T'};
		
		MerTable_t::iterator mi;
		MerTable_t::iterator mjF;
		MerTable_t::iterator mjR;
	
		for (mi = nodes_m.begin(); mi != nodes_m.end(); mi++) {
			Node_t * nodeA = mi->second;
						
			if(nodeA->getTotTmrCov() == 1) { // only process tumor singletons
			
				Mer_t merA = mi->first;
				Mer_t merA_original = merA;
				char old_bp;
				// test changing each bp in the mer
				for (unsigned int i=0; i<merA.size(); i++) {
					
					int qv_cov = (nodeA->cov_distr_tmr[i]).minqv_fwd + (nodeA->cov_distr_tmr[i]).minqv_rev;
					if(qv_cov == 0) { // if loq quality base in tumor
					
						// change bp to any of the 3 other possibile bp
						old_bp = merA[i]; // save old bp
						for (unsigned int j=0; j<4; j++) {
							if(BP[j] != merA[i]) { merA[i] = BP[j]; }
					
							// search for the modified mer (both fwd and rev)
							Mer_t merAr = rc_str(merA);
							mjF = nodes_m.find(merA); // forward
							mjR = nodes_m.find(merAr); // reverse

							// update forward mer (if found)
							if(mjF != nodes_m.end() && mjF != mi) { 
								Mer_t merB = mjF->first;
								Node_t * nodeB = mjF->second;
					
								// merB has only 1 difference from merA
								// only use mers with support >= 2
								if(nodeB->getTotTmrCov() >= MIN_SUPPORT) {
									//cerr << merA_original << "(" << nodeA->getTotTmrCov() << ")\t" << merB << "(" << nodeB->getTotTmrCov() << ")" << "\t" << old_bp << "->" << BP[j] << endl;
						
									// retrive strand info
									unsigned int strand;
									if(nodeA->getTmrCov(FWD)>0) { strand = FWD; }
									else { strand = REV; }
								
									nodeB->incTmrCov(strand);
									nodeB->updateCovDistr((int)(nodeB->getTmrCov(strand)),strand,'T');
									//ref_m->updateCoverage(merB, 'T'); // update referecne k-mer coverage for tumor
									//nodeB->updateCovDistrMinQV(uc_qv,strand,'T');
								}
							}
							
							// update reverse mer (if found)
							if(mjR != nodes_m.end() && mjR != mi) { 
								Mer_t merB = mjR->first;
								Node_t * nodeB = mjR->second;
					
								// merB has only 1 difference from merA
								// only use mers with support >= 2
								if(nodeB->getTotTmrCov() >= MIN_SUPPORT) {
									//cerr << merA_original << "(" << nodeA->getTotTmrCov() << ")\t" << merB << "(" << nodeB->getTotTmrCov() << ")" << "\t" << old_bp << "->" << BP[j] << endl;
						
									// retrive strand info
									unsigned int strand;
									if(nodeA->getTmrCov(FWD)>0) { strand = FWD; }
									else { strand = REV; }
								
									nodeB->incTmrCov(strand);
									nodeB->updateCovDistr((int)(nodeB->getTmrCov(strand)),strand,'T');
									//ref_m->updateCoverage(merB, 'T'); // update referecne k-mer coverage for tumor
									//nodeB->updateCovDistrMinQV(uc_qv,strand,'T');
								}
							}
							
							
						}
						merA[i] = old_bp; //revert changed bp
					}		
				}
			}
		}
	}
	
	bool oneMismatch(Mer_t & a, Mer_t & b) {
		
		int cnt = 0; 
		bool ans = true;
		
		for (unsigned int i=0; i<a.size(); i++) {
			if (a[i]!=b[i]) { cnt++; }
			if(cnt>1) { ans = false; break; }
		}
		
		return ans;
	}
	
};

#endif
