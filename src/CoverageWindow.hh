#ifndef COVERAGEWINDOW_HH
#define COVERAGEWINDOW_HH 1

/****************************************************************************
** Variant.hh
**
** Class for storing coverage information over a genomic window
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

//#include <string>
#include <iostream>
#include <unordered_map>

using namespace std;

class CoverageWindow_t
{
public:

	multiset<int> covset;
	int sum;

	CoverageWindow_t() { sum = 0; }
	
	int getMin() {
		int min = 0;
		if (covset.size()>0) { min = *(covset.begin()); } 
		return min;
	}
	
	int getAvg() { 
		int ans = 0;
		if (covset.size() > 0) {
			ans = (int)floor(sum/covset.size()); 
		}	
		return ans;
	}

	void insert(int c) { covset.insert(c); sum+=c; };
	
	void remove(int c) {
		multiset<int>:: iterator it = covset.find(c);
		if(it != covset.end()) { covset.erase(it); sum-=c; }
	}
	
};

#endif
