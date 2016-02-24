#ifndef ALIGN_HH
#define ALIGN_HH 1

/****************************************************************************
** Align.hh
**
** Gapped alignment based on the Smith-Waterman algorithm
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

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>

void global_align(const std::string & S, const std::string & T,
                  std::string & S_aln, std::string & T_aln,
                  int endfree, int verbose);

void global_align_aff(const std::string & S, const std::string & T,
                      std::string & S_aln, std::string & T_aln,
                      int endfree, int verbose);
 
//void global_cov_align_aff(const std::string & S, const std::string & T, const std::vector<int> & C,
//	                  std::string & S_aln, std::string & T_aln, std::vector<int> & CS_aln, std::vector<int> & CT_aln,
//					  int endfree, int verbose);

#endif
