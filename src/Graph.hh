#ifndef GRAPH_HH
#define GRAPH_HH 1

/****************************************************************************
** Graph.hh
**
** Class for representing and storing a bi-directed de Bruijn graph  
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
#include <deque>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <exception>
#include <stdexcept> // out_of_range exception
#include <iostream>
#include <unordered_map>

#include "align.hh"
#include "util.hh"
#include "Edge.hh"
#include "Node.hh"
#include "Mer.hh"
#include "Ref.hh"
#include "ContigLink.hh"
#include "Path.hh"
#include "Transcript.hh"
#include "Variant.hh"
#include "VariantDB.hh"
#include "CoverageWindow.hh"

#define WHITE 1
#define GREY 2
#define BLACK 3

#define UNDEFINED 0
#define SOURCE 1
#define SINK 2
#define REF 3
#define TMR 4
#define NML 5

using namespace std;

// Graph_t
//////////////////////////////////////////////////////////////////////////

typedef unordered_map<Mer_t, Node_t *> MerTable_t;

class Graph_t
{
public:

	static const char CODE_MAPPED;
	static const char CODE_BASTARD;

	static const string COLOR_ALL;
	static const string COLOR_LOW;
	static const string COLOR_NOVO;
	static const string COLOR_TUMOR;
	static const string COLOR_NORMAL;
	static const string COLOR_SHARED;
	static const string COLOR_SOURCE;
	static const string COLOR_SINK;
	static const string COLOR_TOUCH;

	bool verbose;
	bool VERBOSE;
	bool PRINT_DOT_READS;

	int MIN_QUAL_TRIM;
	int MIN_QUAL_CALL;
	int BUFFER_SIZE;

	int DFS_LIMIT;
	int NODE_STRLEN;

	int K;	
	int MAX_TIP_LEN;
	int MAX_INDEL_LEN;
	int MIN_THREAD_READS;
	int COV_THRESHOLD;
	double MIN_COV_RATIO;
	int LOW_COV_THRESHOLD;
	bool SCAFFOLD_CONTIGS;
	int  INSERT_SIZE;
	int  INSERT_STDEV;
	int  MAX_MISMATCH;
	
	unsigned int NUM_EXTRA_BP; // number of extra bp aftre variant to examine for coverage

	MerTable_t nodes_m;
	int totalreadbp_m;

	Node_t * source_m;
	Node_t * sink_m;

	Ref_t * ref_m;
	bool is_ref_added;

	ReadInfoList_t readid2info;
	int readCycles;
	
	VariantDB_t *vDB; // DB of variants
	Filters filters; // filter thresholds

	Graph_t() : ref_m(NULL), is_ref_added(0), readCycles(0) {
		clear(true); 
		NUM_EXTRA_BP = 20;
	}

	void setDB(VariantDB_t *db) { vDB = db; }
	void setK(int k) { K = k; }
	void setVerbose(bool v) { verbose = v; }
	void setMoreVerbose(bool v) { VERBOSE = v; }
	void setMinQualTrim(int mq) { MIN_QUAL_TRIM = mq; }
	void setMinQualCall(int mq) { MIN_QUAL_CALL = mq; }
	void setBufferSize(int bs) { BUFFER_SIZE = bs; }
	void setDFSLimit(int l) { DFS_LIMIT = l; }
	void setCovThreshold(int ct) { COV_THRESHOLD = ct; }
	void setMinCovRatio(double cr) { MIN_COV_RATIO = cr; }
	void setLowCovThreshold(int ct) { LOW_COV_THRESHOLD = ct; }
	void setPrintDotReads(bool pdr) { PRINT_DOT_READS = pdr; }
	void setNodeStrlen(int ns) { NODE_STRLEN = ns; }
	void setMaxTipLength(int mtl) { MAX_TIP_LEN = mtl; }
	void setMaxIndelLen(int mil) { MAX_INDEL_LEN = mil; }
	void setMinThreadReads(int mtr) { MIN_THREAD_READS = mtr; }
	void setScaffoldContigs(bool sc) { SCAFFOLD_CONTIGS = sc; }
	void setInsertSize(int is) { INSERT_SIZE = is; }
	void setInsertStdev(int is) { INSERT_STDEV = is; }
	void setMaxMismatch(int mm) { MAX_MISMATCH = mm; }
	void setFilters(Filters &fs) { filters = fs; }
	
	bool hasRepeatsInGraphPaths(Ref_t * ref) { return findRepeatsInGraphPaths(source_m, sink_m, F, ref); }

	void clear(bool flag);
	void loadSequence(int readid, const string & seq, const string & qv, bool isRef, int trim5, unsigned int strand);
	void trimAndLoad(int readid, const string & seq, const string & qv, bool isRef, unsigned int strand);
	void trim(int readid, const string & seq, const string & qv, bool isRef);
	void buildgraph(Ref_t * refinfo);

	int countBastardReads();

	int countMappedReads();

	ReadId_t addRead(const string & set, const string & readname, const string & seq, const string & qv, char code, int label, unsigned int strand);

	void addMates(ReadId_t r1, ReadId_t r2);

	void printReads();

	void addPair(const string & set, const string & readname,
		const string & seq1, const string & qv1,
		const string & seq2, const string & qv2,
		char code, int label, unsigned int strand);

	void addUnpaired(const string & set,
		const string & readname,
		const string & seq,
		const string & qv,
		char code,
		int label,
		unsigned int strand);

	void addpaired(const string & set,
		const string & readname,
		const string & seq,
		const string & qv,
		const int mate_id,
		char code,
		int label,
		unsigned int strand);

	//void loadReadsSFA(const string & filename);
	void printAlignment(const string &ref_aln, const string &path_aln, Path_t * path);
	void processPath(Path_t * path, Ref_t * ref, FILE * fp, bool printPathsToFile, int &complete, int &perfect, int &withsnps, int &withindel, int &withmix);
	void processShortPath(Node_t * node, Ref_t * ref, FILE * fp, bool printPathsToFile, int &complete, int &perfect, int &withsnps, int &withindel, int &withmix);
	//Path_t * bfs(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref);
	Path_t * bfs(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref);
	void eka(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref, FILE * fp, bool printPathsToFile);
	void dfs(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref, FILE * fp, bool printPathsToFile);
	bool findRepeatsInGraphPaths(Node_t * source, Node_t * sink, Ori_t dir, Ref_t * ref);
	bool hasTumorOnlyKmer();
	bool hasCycle();
	void hasCycleRec(Node_t * node, Ori_t dir, bool* ans);
	string nodeColor(Node_t * cur, string & who);
	string edgeColor(Node_t * cur, Edge_t & e);
	void printDot(const string & filename, int compid);
	void printFasta(const string & filename);
	void printPairs(const string & filename);
	void markRefEnds(Ref_t * refinfo, int compid);
	void markRefNodes();
	int markConnectedComponents();
	void denovoNodes(const string & filename, const string & refname);
	void alignRefNodes();
	void countRefPath(const string & filename, const string & refname, bool printPathsToFile);
	Node_t * getNode(Mer_t nodeid);
	Node_t * getNode(CanonicalMer_t mer);
	Node_t * getNode(Edge_t & edge);
	void compressNode(Node_t * node, Ori_t dir);
	void compress(int compid);
	void cleanDead();
	void removeNode(Node_t * node);
	void removeLowCov(bool docompression, int compid);
	void removeTips(int compid);
    void greedyTrim();
	void threadReads();
	void checkReadStarts(int compid);
	void updateContigReadStarts();
	void bundleMates();
	void scaffoldContigs();
	void printGraph();
	void printStats(int compid=0);
};

#endif
