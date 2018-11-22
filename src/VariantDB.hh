#ifndef VARIANTDB_HH
#define VARIANTDB_HH 1

/****************************************************************************
** VariantDB.hh
**
** Class for storing multiple variants (DB-style)
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

#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <sstream>
#include "util.hh"
#include "Variant.hh"
#include "sha256.hh"

using namespace std;

struct byPos
{
	bool operator()(pair<string,Variant_t> first, pair<string,Variant_t> second) const {
		
		bool ans = true;
		string chr1 = (first.second).chr; 
		int pos1 = (first.second).pos; 
		string chr2 = (second.second).chr; 
		int pos2 = (second.second).pos; 
	
		int cmp = chr1.compare(chr2);
		if (cmp == 0) { ans = (pos1 < pos2); }
		else if (cmp < 0) { ans = true; }
		else if (cmp > 0) { ans = false; }
	
		return ans;
	}
};

class VariantDB_t
{
public:

	bool LR_MODE;
	map<string,Variant_t> DB; // databbase of variants
	unordered_map<string,int> nCNT; // counts of variants per postion in the normal
	string command_line; // command line used to run the tool

	VariantDB_t(bool lrmode = false) { LR_MODE = lrmode; }
	
	void setLRmode (bool lrmode) { LR_MODE = lrmode; }
	void setCommandLine(string cl) { command_line = cl; }
	void addVar(const Variant_t & v);
	void printHeader(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T);
	void printToVCF(const string version, const string reference, char * date, Filters &fs, string &sample_name_N, string &sample_name_T);
};

#endif
