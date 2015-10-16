#ifndef CONTIGLINK_HH
#define CONTIGLINK_HH 1

/******************************************************************
** ContigLink.hh
**
** Classes for representing and storing links between contigs
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <string>
#include <map>
#include <cmath>

#include "ReadInfo.hh"

using namespace std;

// ContigLink_t
//////////////////////////////////////////////////////////////////////////

class ContigLink_t
{
public:

	ReadId_t rid_m;
	int      linkdist_m;

	ContigLink_t(ReadId_t rid)
		: rid_m(rid), linkdist_m(0)
		{ }

	ContigLink_t(ReadId_t rid, int linkdist)
		: rid_m(rid), linkdist_m(linkdist)
		{ }
};


// ContigLinkList_t
//////////////////////////////////////////////////////////////////////////

class ContigLinkList_t
{
public:

	vector<ContigLink_t> linklist_m;
	unsigned int dupcnt_m;

	ContigLinkList_t()
		: dupcnt_m(0)
		{ }

	void addLink(ReadId_t rid) { linklist_m.push_back(ContigLink_t(rid)); }
	void addLink(ReadId_t rid, int linkdist) { linklist_m.push_back(ContigLink_t(rid, linkdist)); }
	void addDup() { dupcnt_m++; }
	unsigned int dupCnt() { return dupcnt_m; }
	unsigned int linkCnt() { return linklist_m.size(); }

	float mean();
	float stdev(float mean);
};

typedef map<Mer_t, ContigLinkList_t *> ContigLinkMap_t;

#endif

