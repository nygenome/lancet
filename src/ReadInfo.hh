#ifndef READINFO_HH
#define READINFO_HH 1

/****************************************************************************
** ReadInfo.hh
**
** Class for storing general information about a read 
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

#include <string>
#include <vector>

#include "Mer.hh"

#define FWD 1 // forward strand
#define REV 2 // reverse strand

using namespace std;

// ReadId_t
//////////////////////////////////////////////////////////////////////////

typedef int ReadId_t;

// ReadInfo_t
//////////////////////////////////////////////////////////////////////////

class ReadInfo_t
{
public:
	ReadInfo_t(const int label, const string & set, const string & readname, const string & seq, const string & qv, char code, unsigned int strnd, unsigned int mate_order)
		: label_m(label), set_m(set), readname_m(readname), seq_m(seq), qv_m(qv), code_m(code), mateid_m(-1), mate_order_m(mate_order), strand(strnd), trm5(0), trm3(0), isjunk(false)
		{ }

	int            label_m;
	string         set_m;
	string         readname_m;
	string         seq_m;
	string         qv_m;
	char           code_m;
	ReadId_t       mateid_m;
	unsigned short mate_order_m; // is first or second mate? (1=first, 2=mate, 0=unmated)
	Mer_t          contigid_m;
	unsigned int   readstartidx_m;
	unsigned short strand; // FWD or REV
	unsigned short trm5;
	unsigned short trm3;
	bool           isjunk;
};

typedef vector<ReadInfo_t> ReadInfoList_t;

#endif
