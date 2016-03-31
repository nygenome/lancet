#ifndef NODE_HH
#define NODE_HH 1

/***************************************************************************
** Node.hh
**
** Description of a node of the de Bruijn graph 
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
#include <fstream>
#include <string>
#include <assert.h>
#include <unordered_set>

#include "Mer.hh"
#include "Ref.hh"
#include "Edge.hh"
#include "ContigLink.hh"
#include "ReadInfo.hh"
#include "ReadStart.hh"

using namespace std;

// Node_t
//////////////////////////////////////////////////////////////////////////


typedef struct cov_t
{
  int fwd; // total fwd coverage
  int rev; // total rev coverage 
  int minqv_fwd; // min base quality fwd coverage
  int minqv_rev; // min base quality rev coverage
  //int minmq_fwd; // min mapping quality fwd coverage
  //int minmq_rev; // min mapping quality rev coverage
} cov_t;

class Node_t
{
public:

	// members
	//////////////////////////////////////////////////////////////

	Mer_t nodeid_m;
	string qv_m;
	
	string str_m;
	float cov_tmr_m_fwd; // tumor coverage forward
	float cov_tmr_m_rev; // tumor coverage reverse
	float cov_nml_m_fwd; // normal coverage forward
	float cov_nml_m_rev; // normal coverage reverse
	bool isRef_m;
	bool isTumor_m;
	bool isNormal_m;
	bool isSource_m;
	bool isSink_m;
	bool dead_m;
	int  component_m;
	bool touchRef_m;
	int  onRefPath_m;
	int color;
	int MIN_QUAL;

	vector<cov_t> cov_distr_tmr;
	vector<cov_t> cov_distr_nml;
	
	vector<Edge_t> edges_m;
	unordered_set<ReadId_t> reads_m;
	vector<ReadStart_t> readstarts_m;
	ContigLinkMap_t contiglinks_m;
	ReadInfoList_t * readid2info;


	Node_t(Mer_t mer) 
		: nodeid_m(mer), 
		str_m(mer), 
		cov_tmr_m_fwd(0), 
		cov_tmr_m_rev(0), 
		cov_nml_m_fwd(0), 
		cov_nml_m_rev(0),
		isRef_m(false),
		isTumor_m(false),
		isNormal_m(false),
		isSource_m(false),
		isSink_m(false),
		dead_m(false),
		component_m(0),
		touchRef_m(false),
		onRefPath_m(0),
		color(0),
		MIN_QUAL(0)
		{ 	
			cov_distr_tmr.resize(str_m.size()); 
			cov_distr_nml.resize(str_m.size()); 
		}

	friend ostream& operator<<(std::ostream& o, const Node_t & n) { return n.print(o); }
	friend ostream & operator<<(std::ostream & o, const Node_t * n) { return n->print(o); }

	bool isRef() const { return isRef_m; }
	bool isSource() const { return isSource_m; }
	bool isSink() const { return isSink_m; }
	bool isTumor() const { return isTumor_m; }
	bool isNormal() const { return isNormal_m; }
	bool isSpecial() const { return (isSink_m || isSource_m || isRef_m); }
	
	void setIsRef() { isRef_m = true; }
	void setIsSource() { isSource_m = true; }
	void setIsSink() { isSink_m = true; }
	void setIsTumor() { isTumor_m = true; }
	void setIsNormal() { isNormal_m = true; }
	void setMinQV(int q) { MIN_QUAL = q; }
	void setQVstr(string qv) { qv_m = qv; }
	
	void incTmrCov(unsigned int strand) { if(strand == FWD) { cov_tmr_m_fwd++; } if(strand == REV) { cov_tmr_m_rev++; } }
	void incNmlCov(unsigned int strand) { if(strand == FWD) { cov_nml_m_fwd++; } if(strand == REV) { cov_nml_m_rev++; } }
	//void setTmrCov(int c) { cov_tmr_m = c; }
	//void setNmlCov(int c) { cov_nml_m = c; }
	float getTmrCov(unsigned int strand);
	float getNmlCov(unsigned int strand);
	float getTotTmrCov() { return cov_tmr_m_fwd + cov_tmr_m_rev; }
	float getTotNmlCov() { return cov_nml_m_fwd + cov_nml_m_rev; }
	float getTotCov() { return cov_tmr_m_fwd + cov_tmr_m_rev + cov_nml_m_fwd + cov_nml_m_rev; }
	void updateCovDistr(int c, unsigned int strand, char sample);
	void updateCovDistrMinQV(const string & qv, unsigned int strand, char sample);
	void revCovDistr();
	int  minCov();
	int minNon0Cov(char sample);
	int avgCovDistr(char sample);
	
	void setRead2InfoList(ReadInfoList_t * list) { readid2info = list; }
	
	void setColor(int c) { color = c; }
	int getColor() { return color; }

	bool isTandem();
	void addEdge(Mer_t nodeid, Edgedir_t dir, ReadId_t readid);
	void updateEdge(const Mer_t & oldid, Edgedir_t olddir, const Mer_t & newid, Edgedir_t newdir);
	void removeEdge(const Mer_t & nodeid, Edgedir_t dir);
	int getBuddy(Ori_t dir);
	int markRef(Ref_t * ref, int K);
	int degree(Ori_t dir);
	ostream & print(ostream & out) const;
	int strlen() const;
	void addReadStart(ReadId_t readid, int nodeoffset, int trim5, Ori_t ori);
	void revreads();
	void sortReadStarts();
	void addContigLink(Mer_t contigid, ReadId_t rid);
	int cntReadCode(char code);

	int readOverlaps(const Node_t & other);
	
	void swap(cov_t & a, cov_t & b);
};


#endif
