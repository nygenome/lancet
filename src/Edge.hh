#ifndef EDGE_HH
#define EDGE_HH 1

/******************************************************************
** Edge.hh
**
** Class for representing and storing an edge of the de Bruijn graph  
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <set>

#include "Mer.hh"
#include "ReadInfo.hh"


using namespace std;

enum Edgedir_t { FF, FR, RF, RR };

class Edge_t
{
public:

	int flag;
	Mer_t nodeid_m;
	Edgedir_t dir_m;
	vector<ReadId_t> readids_m;


	Edge_t(Mer_t nodeid, Edgedir_t dir)
		: nodeid_m(nodeid), dir_m(dir)
		{ flag = 0; }
	~Edge_t() { };
	
	int getFlag() {return flag; }
	void setFlag(int i) { flag = i; }
	Ori_t startdir() { return edgedir_start(dir_m); }
	Ori_t destdir() { return edgedir_dest(dir_m); }
	string label() { return toString(dir_m) + ":" +  nodeid_m; }

	bool isDir(Ori_t dir);
	int readOverlaps(const Edge_t & other);
	ostream & print(ostream & out) const;
	
	friend ostream& operator<<(std::ostream& o, const Edge_t & e) { return e.print(o); }

    // static functions:

	static char toString(Ori_t dir)
	{
		if (dir == F) { return 'F'; }
		return 'R';
	}

	static Ori_t edgedir_start(Edgedir_t dir)
	{
		if (dir == FF || dir == FR) { return F; }
		return R;
	}

	static Ori_t edgedir_dest(Edgedir_t dir)
	{
		if (dir == FF || dir == RF) { return F; }
		return R;
	}

	static Ori_t flipdir(Ori_t dir)
	{
		if (dir == R) { return F; }
		return R;
	}

	static Edgedir_t flipme(Edgedir_t dir)
	{
		if      (dir == FF) { return RF; }
		else if (dir == FR) { return RR; }
		else if (dir == RF) { return FF; }
		else if (dir == RR) { return FR; }

		return FF;
	}

	static Edgedir_t fliplink(Edgedir_t dir)
	{
		if      (dir == FF) { return RR; }
		else if (dir == FR) { return FR; }
		else if (dir == RF) { return RF; }
		else if (dir == RR) { return FF; }

		return FF;
	}

	static string toString(Edgedir_t dir)
	{
		if      (dir == FF) { return "FF"; }
		else if (dir == FR) { return "FR"; }
		else if (dir == RF) { return "RF"; }
		else if (dir == RR) { return "RR"; }

		return "??";
	}

};

#endif
