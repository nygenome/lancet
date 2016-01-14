#ifndef COVERAGEWINDOW_HH
#define COVERAGEWINDOW_HH 1

/******************************************************************
** Variant.hh
**
** Class for storing coverage infoemation over a genomic window
**
**  Authors: Giuseppe Narzisi
**    Date: January 14, 2016
**
*******************************************************************/

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
