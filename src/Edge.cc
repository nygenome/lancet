#include "Edge.hh"

/******************************************************************
** Edge.cc
**
** Class for representing and storing an edge of the de Bruijn graph  
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

bool Edge_t::isDir(Ori_t dir)
{
	if ((dir == F) && (dir_m == FF || dir_m == FR)) { return true; }
	if ((dir == R) && (dir_m == RR || dir_m == RF)) { return true; }

	return false;
}

int Edge_t::readOverlaps(const Edge_t & other)
{
	set<ReadId_t> myreads;
	for (unsigned int i = 0; i < readids_m.size(); i++)
	{
		myreads.insert(readids_m[i]);
	}

	int retval = 0;
	for (unsigned int j = 0; j < other.readids_m.size(); j++)
	{
		int jj = other.readids_m[j];

		if (myreads.find(jj) != myreads.end())
		{
			retval++;
		}
	}

	return retval;
}

ostream & Edge_t::print(ostream & out) const
{
	out << toString(dir_m) << ":" << nodeid_m;

	out << " [";
	for (unsigned int r = 0; r < readids_m.size(); r++)
	{
		if (r) { out << ","; }
		//out << readid2info[edges_m[i].readids_m[r]].readname_m;
		out << readids_m[r];
	}

	out << "]";

	return out;
}
