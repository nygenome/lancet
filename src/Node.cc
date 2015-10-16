#include "Node.hh"

/******************************************************************
** Node.cc
**
** Implementation of a node of the de Bruijn graph 
** 
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

// isTandem
//////////////////////////////////////////////////////////////

bool Node_t::isTandem()
{
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == nodeid_m)
		{
			return true;
		}
	}

	return false;
}


// addEdge
//////////////////////////////////////////////////////////////

void Node_t::addEdge(Mer_t nodeid, Edgedir_t dir, ReadId_t readid)
{
	if (readid != -1)
	{
		reads_m.insert(readid);
	}

	int edgeid = -1;

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if ((edges_m[i].nodeid_m == nodeid) &&
			(edges_m[i].dir_m == dir))
		{
			if (readid != -1)
			{
				edges_m[i].readids_m.push_back(readid);
			}

			edgeid = i;
			break;
		}
	}

	if (edgeid == -1)
	{
		Edge_t ne(nodeid, dir);

		if (readid != -1)
		{
			ne.readids_m.push_back(readid);
		}
		edges_m.push_back(ne);
	}
}


// updateEdge
//////////////////////////////////////////////////////////////

void Node_t::updateEdge(const Mer_t & oldid, Edgedir_t olddir, 
	const Mer_t & newid, Edgedir_t newdir)
{
	bool found = false;
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == oldid &&
			edges_m[i].dir_m == olddir)
		{
			edges_m[i].nodeid_m = newid;
			edges_m[i].dir_m = newdir;
			return;
		}
	}

	if (!found)
	{
		cerr << "ERROR updating " << nodeid_m << ": didn't find: " << Edge_t::toString(olddir) << ":" << oldid << endl;
		cerr << "wanted to replace with " << Edge_t::toString(newdir) << ":" << newid << endl;
		print(cerr);
	}

	assert(found);
}

// removeEdge
//////////////////////////////////////////////////////////////

void Node_t::removeEdge(const Mer_t & nodeid, Edgedir_t dir)
{
	bool found = false;

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].nodeid_m == nodeid &&
			edges_m[i].dir_m == dir)
		{
			edges_m.erase(edges_m.begin()+i);
			return;
		}
	}

	if (!found)
	{
		cerr << "ERROR removing " << nodeid_m << ": didn't find: " << dir << ":" << nodeid << endl;
	}

	assert(found);
}


// getBuddy: If node has a single edge in the dir direction, return edge id, else -1
//////////////////////////////////////////////////////////////

int Node_t::getBuddy(Ori_t dir)
{
	int retval = -1;

	if (isRef_m) { return retval; }

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].isDir(dir))
		{
			//cerr << edges_m[i] << endl;
			if (retval != -1)
			{
				return -1;
			}

			retval = i;
		}
	}

	// self loop
	if (retval != -1)
	{
		if (edges_m[retval].nodeid_m == nodeid_m)
		{
			return -1;
		}
	}

	return retval;
}

// markRef
//////////////////////////////////////////////////////////////

int Node_t::markRef(Ref_t * ref, int K)
{
	// handle special source/sink nodes
	if (isRef_m) { return 1; }

	CanonicalMer_t cmer;

	touchRef_m = false;

	for (unsigned int i = 0; i < str_m.length()-K+1; i++)
	{
		cmer.set(str_m.substr(i, K));

		if (ref->hasMer(cmer.mer_m))
		{
			touchRef_m = true;
			return 1;
		}
	}

	return 0;
}


// degree
//////////////////////////////////////////////////////////////

int Node_t::degree(Ori_t dir)
{
	int retval = 0;
	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		if (edges_m[i].isDir(dir))
		{
			retval++;
		}
	}

	return retval;
}

// print
//////////////////////////////////////////////////////////////

ostream & Node_t::print(ostream & out) const
{
	out << nodeid_m;

	out << "\t*s\t" << str_m;
	out << "\t*c\t" << cov_m;
	out << "\t*r\t" << isRef_m;

	for (unsigned int i = 0; i < edges_m.size(); i++)
	{
		out << "\t" << edges_m[i];
	}

	return out;
}

// strlen
//////////////////////////////////////////////////////////////

int Node_t::strlen() const
{
	if (isRef_m) { return 0; }
	return str_m.length();
}


// addReadStart
//////////////////////////////////////////////////////////////

void Node_t::addReadStart(ReadId_t readid, int nodeoffset, int trim5, Ori_t ori)
{
	readstarts_m.push_back(ReadStart_t(readid, nodeoffset, trim5, ori));
}

// revreads
//////////////////////////////////////////////////////////////

void Node_t::revreads()
{
	assert(!isRef_m);

	int len = strlen();

	for (unsigned int i = 0; i < readstarts_m.size(); i++)
	{
		ReadStart_t & rs = readstarts_m[i];
		rs.nodeoffset_m = len - 1 - rs.nodeoffset_m;
		rs.ori_m = Edge_t::flipdir(rs.ori_m);
	}
}

// sortReadStarts
//////////////////////////////////////////////////////////////

void Node_t::sortReadStarts()
{
	sort(readstarts_m.begin(), readstarts_m.end(), ReadStart_t::cmpstarts);
}


// addContigLink
//////////////////////////////////////////////////////////////

void Node_t::addContigLink(Mer_t contigid, ReadId_t rid)
{
	ContigLinkMap_t::iterator mi = contiglinks_m.find(contigid);

	if (mi == contiglinks_m.end())
	{
		mi = contiglinks_m.insert(make_pair(contigid, new ContigLinkList_t)).first;
	}

	mi->second->addLink(rid);
}

// cntReadCode
//////////////////////////////////////////////////////////////

int Node_t::cntReadCode(char code)
{
	int retval = 0;

	unordered_set<ReadId_t>::iterator si;
	for (si = reads_m.begin();
	si != reads_m.end();
	si++)
	{
		if ( (*si) < 0 || (*si) > (int)readid2info->size() ) { continue; } // skip over illegal values
		if (readid2info->at(*si).code_m == code)
		{
			retval++;
		}
	}

	return retval;
}

// updateCovDistr
// updated the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::updateCovDistr(int c) 
{
	for (unsigned int i = 0; i < cov_distr.size(); i++)
	{
		cov_distr[i] = c;
	}
}


// revCovDistr
// reverse the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::revCovDistr() 
{
	int tmp;
	int i=0;
	int j=cov_distr.size()-1;
	while(i<j){
		tmp = cov_distr[i];
		cov_distr[i] = cov_distr[j];
		cov_distr[j] = tmp;
		i++;j--;
	}
}

// minCov
// retunr the minimum coverage along the node
//////////////////////////////////////////////////////////////
int Node_t::minCov() 
{
	int min = 10000000;
	for (unsigned int i = 0; i < cov_distr.size(); i++)
	{
		if(cov_distr[i] < min) { min = cov_distr[i]; } 
	}
	return min;
}

// readOverlaps
// return the size of the reads overlap
//////////////////////////////////////////////////////////////
int Node_t::readOverlaps(const Node_t & other)
{
	int retval = 0;
	unordered_set<ReadId_t>::const_iterator it;
	
	for (auto it = other.reads_m.begin(); it != other.reads_m.end(); it++) {

		if (reads_m.find(*it) != reads_m.end())
		{
			retval++;
		}
	}

	return retval;
}
