#include "ContigLink.hh"

/****************************************************************************
** ContigLink.cc
**
** Classes for representing and storing links between contigs
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

float ContigLinkList_t::mean()
{
	float sum = 0;
	float cnt = 0;

	for (unsigned int i = 0; i < linklist_m.size(); i++)
	{
		sum += linklist_m[i].linkdist_m;
		cnt++;
	}

	if (cnt == 0) { return 0.0; }
	return sum / cnt;
}

float ContigLinkList_t::stdev(float mean)
{
	float diffsq = 0;
	float cnt = 0;

	for (unsigned int i = 0; i < linklist_m.size(); i++)
	{
		float diff = linklist_m[i].linkdist_m - mean;
		diffsq += diff*diff;
		cnt++;
	}

	if (cnt == 0) { return 0.0; }
	return sqrtf(diffsq / cnt);
}
