#include "Path.hh"

/****************************************************************************
** Path.cc
**
** Path of a de Bruijn graph 
** Routines to extract sequence and coverage information
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

string Path_t::pathstr()
{
	string retval;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		if (i) { retval += ":"; }

		retval += nodes_m[i]->nodeid_m;
		if (i < edgedir_m.size())
		{
			retval += ":";
			retval += Edge_t::toString(edgedir_m[i]);
		}
	}

	return retval;
}

// pathlen
//////////////////////////////////////////////////////////////

int Path_t::pathlen()
{
	int len = 0;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];

		if (!n->isSpecial())
		{
			len++;
		}
	}

	return len;
}


// str
//////////////////////////////////////////////////////////////

string Path_t::str()
{	
	string retval;
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{				
		Node_t * n = nodes_m[i];
		
		string nstr = n->str_m;

		if (dir == R)
		{
			nstr = rc_str(nstr);
		}

		if (!n->isSpecial())
		{
			if (retval.length() > 0)
			{
				assert(retval.substr(retval.length()-K+1) == nstr.substr(0, K-1));
				retval += nstr.substr(K-1);
			}
			else
			{
				retval = nstr;
			}
		}

		if (i < edgedir_m.size())
		{
			dir = Edge_t::edgedir_dest(edgedir_m[i]);
		}
	}
	
	return retval;
}


// coverage distribution for nodes
//////////////////////////////////////////////////////////////
vector<cov_t> Path_t::covDistr(char sample)
{
	vector<cov_t> path_coverage;
	vector<cov_t> node_coverage;
	vector<cov_t> C;
	
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);
	
	//cerr << "Num nodes in path: " << nodes_m.size() << endl;
	
	path_coverage.clear();	
	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		node_coverage.clear();
		Node_t * n = nodes_m[i]; 
				
		if(sample == 'T') { C = n->cov_distr_tmr; }
		else if(sample == 'N') { C = n->cov_distr_nml; }
		else { cerr << "Error: unrecognized sample " << sample << endl; }
				
		if (dir == R) {
			for (unsigned int j=C.size(); j>0; j--) { 
				node_coverage.push_back(C[j-1]); 
			}
		}
		else {
			for (unsigned int j=0; j < C.size(); j++) { 
				node_coverage.push_back(C[j]); 
			}
		}
		
		if (!n->isSpecial()) {	
						
			if (path_coverage.size() == 0) { // first node
								
				for (unsigned int j=0; j < node_coverage.size(); j++) { 
					path_coverage.push_back(node_coverage[j]); 
				}
			}
			else { // not the first node: update coverage of overlapping region
				// add coverage info for the new base-pairs 
				
				/*
				int p = (path_coverage.size())-K+1;
				assert(p>0); 
				for (int l = 0; l < K-1; l++) {
					// tumor
					int path_cov = path_coverage[p+l].fwd + path_coverage[p+l].rev;
					int node_cov = node_coverage[l].fwd + node_coverage[l].rev;
					if(path_cov < node_cov) { path_coverage[p+l] = node_coverage[l]; }					
				}
				*/
				for (unsigned int j = (K-1); j < node_coverage.size(); j++) {
					path_coverage.push_back(node_coverage[j]);
				}				
			}
		}
		
		if (i < edgedir_m.size())
		{
			dir = Edge_t::edgedir_dest(edgedir_m[i]);
		}
	}

	return path_coverage;
}

// coverage at position
//////////////////////////////////////////////////////////////
/*
int Path_t::covAt(int pos, char sample, unsigned int strand)
{
	int retval = -1; 
	int p = 0;
	vector<int> C;
	vector<int> coverage;
	
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);
	
	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		coverage.clear();
		Node_t * n = nodes_m[i];
		
		if(sample == 'T') {  // tumor coverage
			if(strand == FWD) { C = n->cov_distr_tmr_fwd; } 
			if(strand == REV) { C = n->cov_distr_tmr_rev; } 			
		}
		else if(sample == 'N') { // normal coverage
			if(strand == FWD) { C = n->cov_distr_nml_fwd; } 
			if(strand == REV) { C = n->cov_distr_nml_rev; } 	
		}  
		
		if (dir == R) {
			for (unsigned int j=C.size(); j>0; j--) { 
				coverage.push_back(C[j-1]);
			}
		}
		else {
			for (unsigned int j=0; j < C.size(); j++) { 
				coverage.push_back(C[j]); 
			}
		}
		
		if (!n->isSpecial())
		{
			unsigned int j = 0;
			if (p > 0) { // if not first node, scan only the extra base-pairs for 
				j = K-1;
			}
			for (; j < coverage.size(); j++) {
				//if(p == pos) { return n->cov_distr[j]; }
				if(p == pos) { return coverage[j]; }
				p++;
			}
		}
		
		if (i < edgedir_m.size())
		{
			dir = Edge_t::edgedir_dest(edgedir_m[i]);
		}
	}

	return retval;
}
*/

// coverage distribution for edges
//////////////////////////////////////////////////////////////

vector<float> Path_t::readCovNodes()
{
	vector<float> nodes_coverage;
		
	//cerr << "Num Nodes in path: " << nodes_m.size() << endl;
	
	nodes_coverage.clear();	
	for (unsigned int i = 1; i < nodes_m.size(); i++)
	{
		Node_t* n1 = nodes_m[i-1];
		Node_t* n2 = nodes_m[i]; 
		
		nodes_coverage.push_back(n1->readOverlaps(*n2)); 
	}

	return nodes_coverage;
}

// cov
//////////////////////////////////////////////////////////////

float Path_t::cov(char sample)
{
	float covsum = 0;
	float strlen = 0;
	float C = 0;
	
	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];
		
		if(sample == 'T') { C = n->getTotTmrCov(); } // tumor coverage
		else if(sample == 'N') { C = n->getTotNmlCov(); } // normal coverage
		else if(sample == 'A') { C = n->getTotTmrCov() + n->getTotNmlCov(); } // normal + tumor coverage

		if (!n->isSpecial())
		{
			int merlen = n->strlen() - K + 1;
			covsum += C * merlen;
			strlen += merlen;
		}
	}

	return covsum / strlen;
}

// mincov
//////////////////////////////////////////////////////////////

float Path_t::mincov(char sample)
{
	float mincov = -1;
	float C = 0;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];
		
		if(sample == 'T') { C = n->getTotTmrCov(); } // tumor coverage
		else if(sample == 'N') { C = n->getTotNmlCov(); } // normal coverage
		else if(sample == 'A') { C = n->getTotTmrCov() + n->getTotNmlCov(); } // normal + tumor coverage

		if (!n->isSpecial())
		{
			if ((mincov == -1) || (C < mincov))
			{
				mincov = C;
			}
		}
	}

	return mincov;
}

// maxcov
//////////////////////////////////////////////////////////////

float Path_t::maxcov(char sample)
{
	float maxcov = -1;
	float C = 0;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];
		
		if(sample == 'T') { C = n->getTotTmrCov(); } // tumor coverage
		else if(sample == 'N') { C = n->getTotNmlCov(); } // normal coverage
		else if(sample == 'A') { C = n->getTotTmrCov() + n->getTotNmlCov(); } // normal + tumor coverage

		if (!n->isSpecial())
		{
			if ((maxcov == -1) || (C > maxcov))
			{
				maxcov = C;
			}
		}
	}

	return maxcov;
}

// pathcontig
//////////////////////////////////////////////////////////////

Node_t * Path_t::pathcontig(int pos)
{
	int curpos = 0;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];

		if (!n->isSpecial())
		{
			int span = n->str_m.length();

			if (curpos + span >= pos)
			{
				// in the right node
				return n;
			}

			curpos += span - K + 1;
		}
	}

	return NULL;
}

// contains
//////////////////////////////////////////////////////////////

int Path_t::hasCycle(Node_t * node)
{
	if (hasCycle_m)
		return hasCycle_m;

	for (vector<Node_t *>::iterator ni = nodes_m.begin(); ni != nodes_m.end(); ni++)
	{
		if (*ni == node)
		{
			hasCycle_m = 1;
			return 1;
		}
	}

	return 0;
}

// return true if tumor only node is found in the path
//////////////////////////////////////////////////////////////

bool Path_t::hasTumorOnlyNode()
{
	bool ans = false;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];
		
		if (n->isTumor() && !n->isNormal()) { 
			ans = true; 
			break; // exit as soon as tumor spcific node is found
		}
	}

	return ans;
}
