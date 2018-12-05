#include "Node.hh"

/****************************************************************************
** Node.cc
**
** Implementation of a node of the de Bruijn graph
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


// addBX
// add 10x barcode to set of barcodes for this node
// return false if the insertion was not succesfull
//////////////////////////////////////////////////////////////
bool Node_t::addBX(std::string & bx, unsigned int strand, int label) { 

	pair<unordered_set<string>::iterator,bool> it;
	if (label == TMR) {
		if(strand == FWD) { it = bxset_tmr_fwd.insert(bx); }
		if(strand == REV) { it = bxset_tmr_rev.insert(bx); }
	}
	if (label == NML) {
		if(strand == FWD) { it = bxset_nml_fwd.insert(bx); }
		if(strand == REV) { it = bxset_nml_rev.insert(bx); }
	}
	return it.second;
}

// addHP
// add 10x haplotype number to set of haplotypes for this node
// return false if the insertion was not succesfull
//////////////////////////////////////////////////////////////
void Node_t::addHP(int hp, int label) { 
	if (label == TMR) { hpset_tmr.at(hp) += 1; }
	if (label == NML) { hpset_nml.at(hp) += 1; }
}

// hasBX
// return true if the barcode is already present in the bxset
//////////////////////////////////////////////////////////////
bool Node_t::hasBX(std::string & bx, int label) {

	bool ans = false;
	
	if(label == TMR) {
		auto got_fwd = bxset_tmr_fwd.find(bx);
		if ( got_fwd != bxset_tmr_fwd.end() ) { ans = true; }
		else {
			auto got_rev = bxset_tmr_rev.find(bx);
			if ( got_rev != bxset_tmr_rev.end() ) { ans = true; }
		}
	}
	
	if(label == NML) {
		auto got_fwd = bxset_nml_fwd.find(bx);
		if ( got_fwd != bxset_nml_fwd.end() ) { ans = true; }
		else {
			auto got_rev = bxset_nml_rev.find(bx);
			if ( got_rev != bxset_nml_rev.end() ) { ans = true; }
		}
	}
	
	return ans;
}

// BXcnt
// return BX coverage by strand
//////////////////////////////////////////////////////////////
int Node_t::BXcnt(unsigned int strand, int label) {
	
	int cnt = -1;
	
	if (label == TMR) {
		if(strand == FWD) { cnt = bxset_tmr_fwd.size(); }
		if(strand == REV) { cnt = bxset_tmr_rev.size(); }
	}
	
	if (label == NML) {
		if(strand == FWD) { cnt = bxset_nml_fwd.size(); }
		if(strand == REV) { cnt = bxset_nml_rev.size(); }
	}
	
	return cnt;
}

// HPcnt
// return haplotype coverage
//////////////////////////////////////////////////////////////
int Node_t::HPcnt(unsigned int hp_num, int label) {
	
	int cnt = -1;
	
	if (label == TMR) { cnt = hpset_tmr.at(hp_num); }
	if (label == NML) { cnt = hpset_nml.at(hp_num); }
	
	return cnt;
}

// isTandem
//////////////////////////////////////////////////////////////

bool Node_t::isTandem()
{
	for (unsigned int i = 0; i < edges_m.size(); ++i)
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

	for (unsigned int i = 0; i < edges_m.size(); ++i)
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
	for (unsigned int i = 0; i < edges_m.size(); ++i)
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

	for (unsigned int i = 0; i < edges_m.size(); ++i)
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

	//if (isRef()) { return retval; }
	if (isSpecial()) { return retval; }

	for (unsigned int i = 0; i < edges_m.size(); ++i)
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
	if (isSource() || isSink()) { return 1; }

	CanonicalMer_t cmer;

	touchRef_m = false;

	for (unsigned int i = 0; i < str_m.length()-K+1; ++i)
	{
		Mer_t m = str_m.substr(i, K);
		cmer.set(m);
		//cmer.set(str_m.substr(i, K));

		if (ref->hasMer(cmer.mer_m))
		{
			touchRef_m = true;
			//isRef_m = true;
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
	for (unsigned int i = 0; i < edges_m.size(); ++i)
	{
		if (edges_m[i].isDir(dir))
		{
			++retval;
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
	out << "\t*c\t" << cov_tmr_m_fwd;
	out << "\t*c\t" << cov_nml_m_fwd;
	out << "\t*c\t" << cov_tmr_m_rev;
	out << "\t*c\t" << cov_nml_m_rev;
	out << "\t*r\t" << isRef();

	for (unsigned int i = 0; i < edges_m.size(); ++i)
	{
		out << "\t" << edges_m[i];
	}

	return out;
}

// strlen
//////////////////////////////////////////////////////////////

int Node_t::strlen() const
{
	//if (isRef()) { return 0; }
	if (isSpecial()) { return 0; }
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
	//assert(!isRef());
	assert(!isSpecial());

	int len = strlen();

	for (unsigned int i = 0; i < readstarts_m.size(); ++i)
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

	for (auto si = reads_m.begin();
	si != reads_m.end();
	++si)
	{
		if ( (*si) < 0 || (*si) > (int)readid2info->size() ) { continue; } // skip over illegal values
		if (readid2info->at(*si).code_m == code)
		{
			++retval;
		}
	}

	return retval;
}

// isstatusCnt
// return true if mor than 90% of the postion in the nodes are 
// of type c (T or N)
//////////////////////////////////////////////////////////////
bool Node_t::isStatusCnt(char c) {
	
	bool ans = false;
	int cnt = 0;
	unsigned int N = 0;
	
	for (unsigned int i = (K-1); i < cov_status.size(); ++i) {
		//cerr << cov_status[i] << " ";
		++N;
		if(cov_status[i] == c) { ++cnt; }
	}
	
	double prcnt = (double(cnt)/double(N));
	if ( prcnt > 0.8) { ans = true; }
	//cerr << "Percent tumor: " << prcnt << " " << cnt << "/" << N << endl;
	
	return ans;
}

// updateCovStatus
// updated the coverage status along the node string
//////////////////////////////////////////////////////////////
void Node_t::updateCovStatus(char c) 
{
	for (unsigned int i = 0; i < cov_status.size(); ++i) {
		
		/*
		switch(cov_status[i]) {
		    case 'E': cov_status[i] = c;
					break;
		    case 'T': if(c == 'N') { cov_status[i] = 'B'; }
					break;
		 	case 'N': if(c == 'T') { cov_status[i] = 'B'; }
					break;
			default : cov_status[i] = c;
					break;
		}
		*/
		if(cov_status[i] == 'E') { cov_status[i] = c; }
		else if(cov_status[i] != c) { cov_status[i] = 'B'; }
		else { cov_status[i] = c; }
	}
}

// updateCovDistr
// updated the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::updateCovDistr(int cov, const string & qv, unsigned int strand, int sample) 
{
	vector<cov_t> * cov_distr = NULL;
	
	if(sample == TMR)      { cov_distr = &cov_distr_tmr; }
	else if(sample == NML) { cov_distr = &cov_distr_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
 	//string::const_iterator it=qv.begin();
	for (unsigned int i = 0; i < cov_distr->size(); ++i) {
		
		if(strand == FWD) { 
			((*cov_distr)[i]).fwd = cov;
			if(qv[i] >= MIN_QUAL) { 
				++(((*cov_distr)[i]).minqv_fwd);
			}
		}
		else if(strand == REV) { 
			((*cov_distr)[i]).rev = cov;
			if(qv[i] >= MIN_QUAL) { 
				++(((*cov_distr)[i]).minqv_rev);
			}
		}
		//if (it!=qv.end()) { it++; }
		//else {cerr << "Error: reached end of quality string (qv)" << endl; }
	}
}

// updateHPCovDistr
// updated the haplotype coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::updateHPCovDistr(int hp0_cov, int hp1_cov, int hp2_cov, int sample) 
{
	vector<cov_t> * cov_distr = NULL;
	
	if(sample == TMR)      { cov_distr = &cov_distr_tmr; }
	else if(sample == NML) { cov_distr = &cov_distr_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	for (unsigned int i = 0; i < cov_distr->size(); ++i) {
		((*cov_distr)[i]).hp0 = hp0_cov;
		((*cov_distr)[i]).hp1 = hp1_cov;
		((*cov_distr)[i]).hp2 = hp2_cov;
	}
}

// avgCovDistr
// average coverage of non-zero elements
//////////////////////////////////////////////////////////////
int Node_t::avgCovDistr(char sample)
{
	vector<cov_t> cov_distr;
	
	if(sample == 'T')      { cov_distr = cov_distr_tmr; }
	else if(sample == 'N') { cov_distr = cov_distr_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	int sum = 0;
	int cnt = 0;
	for (unsigned int i = 0; i < cov_distr.size(); ++i)
	{
		int totcov = cov_distr[i].fwd + cov_distr[i].rev;
		if(totcov !=0) {
			sum += totcov;
			++cnt;
		}
	}
	//cerr << "SPANNER: (" << sum << "," << cnt << ")" << endl; 

	int ans = 0;
	if (sum>0) { ans = floor(float(sum)/float(cnt)); }
	
	return ans;
}

// swap two cov_t
//////////////////////////////////////////////////////////////////////////

void Node_t::swap(cov_t & a, cov_t & b) 
{
	cov_t tmp = a;
	a = b;
	b = tmp;
}

// revCovDistr
// reverse the coverage distribution along the node string
//////////////////////////////////////////////////////////////
void Node_t::revCovDistr() 
{
	int i=0;
	int j=cov_distr_tmr.size()-1;
	while(i<j){
		swap(cov_distr_tmr[i], cov_distr_tmr[j]);
		swap(cov_distr_nml[i], cov_distr_nml[j]);
		++i;--j;
	}
}

// minCov
// retunr the minimum coverage along the node
//////////////////////////////////////////////////////////////
int Node_t::minNon0Cov(char sample) 
{
	
	vector<cov_t> cov_distr;

	if(sample == 'T')      { cov_distr = cov_distr_tmr; }
	else if(sample == 'N') { cov_distr = cov_distr_nml; }
	else { cerr << "Error: unrecognized sample " << sample << endl; }
	
	int min = 10000000;
	for (unsigned int i = 0; i < cov_distr.size(); ++i)
	{
		int totcov = cov_distr[i].fwd + cov_distr[i].rev;
		if( (totcov > 0) && (totcov < min) ) { min = totcov; }
	}
	
	return min;
}

// computeMinCov
// return the minimum coverage along the node
//////////////////////////////////////////////////////////////
void Node_t::computeMinCov()
{
	int min = 10000000;
	int minQV = 10000000;
	for (unsigned int i = 0; i < cov_distr_tmr.size(); ++i)
	{
		int totcov = cov_distr_tmr[i].fwd + cov_distr_tmr[i].rev + cov_distr_nml[i].fwd + cov_distr_nml[i].rev;
		int totcovQV = cov_distr_tmr[i].minqv_fwd + cov_distr_tmr[i].minqv_rev + cov_distr_nml[i].minqv_fwd + cov_distr_nml[i].minqv_rev;
		
		if(totcov < min) { min = totcov; } 
		if(totcovQV < minQV) { minQV = totcovQV; } 
		
	}
	mincov = min;
	mincovQV = minQV;
}

// readOverlaps
// return the size of the reads overlap
//////////////////////////////////////////////////////////////
int Node_t::readOverlaps(const Node_t & other)
{
	int retval = 0;
	
	for (auto it = other.reads_m.begin(); it != other.reads_m.end(); ++it) {

		if (reads_m.find(*it) != reads_m.end())
		{
			++retval;
		}
	}

	return retval;
}

// hasOverlappingMate
// return true if the k-mer comes from the same fragment (overlapping mates)
//////////////////////////////////////////////////////////////
bool Node_t::hasOverlappingMate(string & read_name, int id)
{	
	bool ans = false;
	
	if(id == 1) {
		//if (mate2_name.find(read_name) != mate2_name.end()) { ans = true; }
	    if (binary_search (mate2_name.begin(), mate2_name.end(), read_name)) { ans = true; }		

		//for (vector<string>::iterator it2 = mate2_name.begin() ; it2 != mate2_name.end(); ++it2) {
		//	if ((*it2) == read_name) { ans = true; }
		//}
	}
	
	if(id == 2) {
		//if (mate1_name.find(read_name) != mate1_name.end()) { ans = true; }
	    if (binary_search (mate1_name.begin(), mate1_name.end(), read_name)) { ans = true; }		
			
	    //for (vector<string>::iterator it1 = mate1_name.begin() ; it1 != mate1_name.end(); ++it1) {
		//	if ((*it1) == read_name) { ans = true; }
		//}
	}
	
	return ans;
}

// add mate name to the set of mates containing this kmer
// also store  mate order (1st or 2nd in pair) 
void Node_t::addMateName(string & read_name, int id) 
{	
	//if(id == 1) { mate1_name.insert(read_name); }
	//if(id == 2) { mate2_name.insert(read_name); }
	if(id == 1) { mate1_name.push_back(read_name); }
	if(id == 2) { mate2_name.push_back(read_name); }
}


// return  coverage on the input strand for tumor or normal sample
//////////////////////////////////////////////////////////////
float Node_t::getCov(unsigned int strand, int label) { 
	float ans = 0;
	
	if (label == TMR) {
		if(strand == FWD) { ans = cov_tmr_m_fwd; } 
		if(strand == REV) { ans = cov_tmr_m_rev; } 
	}	
	if (label == NML) {
		if(strand == FWD) { ans = cov_nml_m_fwd; } 
		if(strand == REV) { ans = cov_nml_m_rev; } 
	}
	
	return ans;
}

// increase node coverage for tumor or normal
//////////////////////////////////////////////////////////////
void Node_t::incCov(unsigned int strand, int label) { 
	
	if (label == TMR) {
		if(strand == FWD) { cov_tmr_m_fwd++; } 
		else if(strand == REV) { cov_tmr_m_rev++; }
	}
	if (label == NML) {
		if(strand == FWD) { cov_nml_m_fwd++; } 
		else if(strand == REV) { cov_nml_m_rev++; }
	}
}
