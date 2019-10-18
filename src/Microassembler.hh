#ifndef MICROASSEMBLER_HH
#define MICROASSEMBLER_HH 1

/****************************************************************************
** Microassembler.hh
**
** Tool for localized assembly of genomic regions using
** the de Bruijn paradigm to detect genetic variants
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

#include <unistd.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <assert.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <getopt.h>
#include <pthread.h>
#include <time.h>       /* time_t, struct tm, time, localtime, strftime */

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/SamReadGroupDictionary.h"
#include "api/SamReadGroup.h"
#include "htslib/faidx.h"

#include "align.hh"
#include "util.hh"
#include "Mer.hh"
#include "Ref.hh"
#include "ReadInfo.hh"
#include "ReadStart.hh"
#include "Transcript.hh"
#include "Edge.hh"
#include "ContigLink.hh"
#include "Node.hh"
#include "Path.hh"
#include "Graph.hh"
#include "VariantDB.hh"
#include "ErrorCorrector.hh"

using namespace std;
using namespace HASHMAP;
using namespace BamTools;

#define bvalue(value) ((value ? "true" : "false"))

class Microassembler {

public:
	
	int ID;
	int BUFFER_SIZE;
	int WINDOW_SIZE;
	
	// Configuration
	//////////////////////////////////////////////////////////////////////////

	bool verbose;
	bool VERBOSE;
	bool PRINT_DOT_READS;
	bool KMER_RECOVERY;
	bool PRINT_ALL;

	int MIN_QV_CALL;
	int MIN_QV_TRIM;
	int QV_RANGE;
	int MIN_QUAL_CALL;
	int MIN_QUAL_TRIM;
	int MIN_MAP_QUAL;
	int MAX_DELTA_AS_XS;

	string TUMOR;
	string NORMAL;
	string RG_FILE;
	string REFFILE;
	string READSET;

	//string PREFIX;

	int minK;
	int maxK;
	int MAX_TIP_LEN;
	unsigned int MIN_THREAD_READS;
	int COV_THRESHOLD;
	double MIN_COV_RATIO;
	int LOW_COV_THRESHOLD;
	int MAX_AVG_COV;
	
	//STR parameters
	int MAX_UNIT_LEN;
	int MIN_REPORT_UNITS;
	int MIN_REPORT_LEN;
	int DIST_FROM_STR;
	
	bool LR_MODE;
	bool PRIMARY_ALIGNMENT_ONLY;
	bool XA_FILTER;
	bool ACTIVE_REGION_MODULE;
	bool SCAFFOLD_CONTIGS;
	int  INSERT_SIZE;
	int  INSERT_STDEV;

	int  NODE_STRLEN;

	int DFS_LIMIT;
	int MAX_INDEL_LEN;
	int MAX_MISMATCH;
		
	Filters * filters; // filter thresholds
	
	string sample_name_normal;
	string sample_name_tumor;
	
	// data structures
	//////////////////////////////////////////////////////////////////////////
	
	int graphCnt;
	int num_skip; // number of regions/windows skipped 
	set<string> readgroups;
	set<string> RG_father;
	set<string> RG_mother;
	set<string> RG_self;
	set<string> RG_sibling;
	
	map<string, Ref_t *> * reftable; // table of references to analyze
	VariantDB_t vDB; // variants DB
	
	int num_snv_only_regions;
	int num_indel_only_regions;
	int num_softclip_only_regions;
	int num_indel_or_softclip_regions;
	int num_snv_or_indel_regions;
	int num_snv_or_softclip_regions;
	int num_snv_or_indel_or_softclip_regions;
	
	Microassembler(bool lrmode) { 
		
		LR_MODE = lrmode;
		vDB.setLRmode(lrmode);
		
		graphCnt = 0;
		num_skip = 0;
		
		ACTIVE_REGION_MODULE = true;
		PRIMARY_ALIGNMENT_ONLY = false;
		XA_FILTER = false;
		
		BUFFER_SIZE = 10*1024;
		WINDOW_SIZE = 600;

		verbose			= false;
		VERBOSE         = false;
		PRINT_DOT_READS = true;
		KMER_RECOVERY	= false;
		PRINT_ALL       = false;

		MIN_QV_CALL    = 10;
		MIN_QV_TRIM    = 10;
		QV_RANGE       = '!';
		MIN_QUAL_CALL  = MIN_QV_CALL + QV_RANGE;
		MIN_QUAL_CALL  = MIN_QV_TRIM + QV_RANGE;
		MIN_MAP_QUAL   = 0;
		MAX_DELTA_AS_XS = 5;

		READSET = "qry";

		minK = 11;
		maxK = 101;
		MAX_TIP_LEN = minK;
		MIN_THREAD_READS = 3;
		COV_THRESHOLD = 5;
		MIN_COV_RATIO = 0.01;
		LOW_COV_THRESHOLD = 1;
		MAX_AVG_COV = 10000;

		SCAFFOLD_CONTIGS = 0;
		INSERT_SIZE = 150;
		INSERT_STDEV = 15;

		NODE_STRLEN = 100;

		RG_FILE = "";

		DFS_LIMIT = 1000000;
		MAX_INDEL_LEN = 500;
		MAX_MISMATCH = 2;
		
		num_snv_only_regions = 0;
		num_indel_only_regions = 0;
		num_softclip_only_regions = 0;
		num_indel_or_softclip_regions = 0;
		num_snv_or_indel_regions = 0;
		num_snv_or_softclip_regions = 0;
		num_snv_or_indel_or_softclip_regions = 0;
	}
		
	~Microassembler() { }
	
	void loadRefs(const string & filename);
	void loadRG(const string & filename, int member);
	int processGraph(Graph_t & g, const string & refname, int minK, int maxK);
	int run(int argc, char** argv);
	bool extractReads(BamReader &reader, Graph_t &g, Ref_t *refinfo, BamRegion &region, int &readcnt, int code);
	bool isActiveRegion(BamReader &reader, Ref_t *refinfo, BamRegion &region, int code);
	int processReads();
	void setFilters(Filters * fs) { filters = fs; vDB.setFilters(fs); }
	void setID(int i) { ID = i; }
	string retriveSampleName(SamHeader &header);
};

#endif
