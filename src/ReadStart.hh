#ifndef READSTART_HH
#define READSTART_HH 1

/****************************************************************************
** ReadStart.hh
**
** Class for storing start location of a read in sequence path
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

class ReadStart_t
{
public:
	ReadStart_t(ReadId_t readid, int nodeoffset, int trim5, Ori_t ori)
		: readid_m(readid), nodeoffset_m(nodeoffset), trim5_m(trim5), ori_m(ori)
		{ }

	static bool cmpstarts (const ReadStart_t & a, const ReadStart_t & b) 
		{ return a.nodeoffset_m < b.nodeoffset_m; }

	ReadId_t readid_m;
	int      nodeoffset_m;
	int      trim5_m;
	Ori_t    ori_m;
};

#endif
