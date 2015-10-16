#include "Path.hh"

/******************************************************************
** Path.cc
**
** Path of a de Bruijn graph 
** Routines to extract sequence and coverage information
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

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

		if (!n->isRef_m)
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

		if (!n->isRef_m)
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

// coverage distribution for nodes in string format
//////////////////////////////////////////////////////////////

string Path_t::covstr() {
	stringstream ss;
	vector<float> node_coverage;
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		node_coverage.clear();
		Node_t * n = nodes_m[i]; 
				
		if (dir == R) {
			for (unsigned int j=n->cov_distr.size(); j>0; j--) { 
				node_coverage.push_back(n->cov_distr[j-1]); 
			}
		}
		else {
			for (unsigned int j=0; j < n->cov_distr.size(); j++) { 
				node_coverage.push_back(n->cov_distr[j]); 
			}
		}
		
		if (!n->isRef_m) {		
		
			if ((ss.str()).size() == 0) { // first node
				for (unsigned int j=0; j < node_coverage.size(); j++) { 
					if (i != (nodes_m.size()-1)) { ss << node_coverage[j] << " "; }
					else { ss << node_coverage[j]; }
				}
			}
			else { // not the first node: update coverage of overlapping region
				// add coverage info for the new base-pairs 	
				for (unsigned int j = (K-1); j < node_coverage.size(); j++) {
					if (i != (nodes_m.size()-1)) { ss << node_coverage[j] << " "; }
					else { ss << node_coverage[j]; }
				}
			}
		}
		
		if (i < edgedir_m.size())
		{
			dir = Edge_t::edgedir_dest(edgedir_m[i]);
		}
	}

	return ss.str();
}

// coverage distribution for nodes
//////////////////////////////////////////////////////////////

vector<int> Path_t::covDistr()
{
	vector<int> path_coverage;
	vector<int> node_coverage;
	
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);
	
	//cerr << "Num nodes in path: " << nodes_m.size() << endl;

	path_coverage.clear();	
	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		node_coverage.clear();
		Node_t * n = nodes_m[i]; 
		
		if (dir == R) {
			for (unsigned int j=n->cov_distr.size(); j>0; j--) { 
				node_coverage.push_back(n->cov_distr[j-1]); 
			}
		}
		else {
			for (unsigned int j=0; j < n->cov_distr.size(); j++) { 
				node_coverage.push_back(n->cov_distr[j]); 
			}
		}
		
		if (!n->isRef_m) {		
		
			if (path_coverage.size() == 0) { // first node
				for (unsigned int j=0; j < node_coverage.size(); j++) { 
					path_coverage.push_back(node_coverage[j]); 
				}
			}
			else { // not the first node: update coverage of overlapping region
				// add coverage info for the new base-pairs 	
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

int Path_t::covAt(int pos)
{
	int retval = -1; 
	int p = 0;
	vector<int> coverage;
	
	Ori_t dir = Edge_t::edgedir_start(edgedir_m[0]);
	
	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		coverage.clear();
		Node_t * n = nodes_m[i];
		
		if (dir == R) {
			for (unsigned int j=n->cov_distr.size(); j>0; j--) { 
				coverage.push_back(n->cov_distr[j-1]);
			}
		}
		else {
			for (unsigned int j=0; j < n->cov_distr.size(); j++) { 
				coverage.push_back(n->cov_distr[j]); 
			}
		}
		
		if (!n->isRef_m)
		{
			unsigned int j = 0;
			if (p > 0) { // if not first node, scan only the extra base-pairs for 
				j = K-1;
			}
			//for (; j < n->cov_distr.size(); j++) {
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

float Path_t::cov()
{
	float covsum = 0;
	float strlen = 0;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];

		if (!n->isRef_m)
		{
			int merlen = n->strlen() - K + 1;
			covsum += n->cov_m * merlen;
			strlen += merlen;
		}
	}

	return covsum / strlen;
}

// mincov
//////////////////////////////////////////////////////////////

float Path_t::mincov()
{
	float mincov = -1;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];

		if (!n->isRef_m)
		{
			if ((mincov == -1) || (n->cov_m < mincov))
			{
				mincov = n->cov_m;
			}
		}
	}

	return mincov;
}

// maxcov
//////////////////////////////////////////////////////////////

float Path_t::maxcov()
{
	float maxcov = -1;

	for (unsigned int i = 0; i < nodes_m.size(); i++)
	{
		Node_t * n = nodes_m[i];

		if (!n->isRef_m)
		{
			if ((maxcov == -1) || (n->cov_m > maxcov))
			{
				maxcov = n->cov_m;
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

		if (!n->isRef_m)
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
