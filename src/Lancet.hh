#ifndef LANCET_HH
#define LANCET_HH 1

/****************************************************************************
** Lancet.hh
**
** Tool for localized colored de Bruijn graph assembly to detect 
** somatic genetic variants
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

#include "Microassembler.hh"

string VERSION = "1.1.0, October 18 2019";
string COMMAND_LINE;

/****  configuration parameters ****/
int NUM_THREADS = 1;

bool LR_MODE = false; // linked-reads mode
bool XA_FILTER = false;
bool PRIMARY_ALIGNMENT_ONLY = false;
bool ACTIVE_REGIONS = true;
bool verbose = false;
bool VERBOSE = false;
bool KMER_RECOVERY = false;
bool PRINT_ALL = false;
bool PRINT_DOT_READS = true;
int MIN_QV_TRIM = 10;
int MIN_QV_CALL = 17;
int QV_RANGE = '!';
int MIN_QUAL_TRIM = MIN_QV_TRIM + QV_RANGE;
int MIN_QUAL_CALL = MIN_QV_CALL + QV_RANGE;
int MIN_MAP_QUAL = 15;
int MAX_DELTA_AS_XS = 5;
int WINDOW_SIZE = 600;
int PADDING = 250;

string TUMOR;
string NORMAL;
string RG_FILE;
string REFFILE;
string BEDFILE;
string REGION;

int minK = 11;
int maxK = 101;
int MAX_TIP_LEN = minK;
unsigned int MIN_THREAD_READS = 1;
int COV_THRESHOLD = 5;
double MIN_COV_RATIO = 0.01;
int LOW_COV_THRESHOLD = 1;
int MAX_AVG_COV = 10000;
int NODE_STRLEN = 100;
int DFS_LIMIT = 1000000;
int MAX_INDEL_LEN = 500;
int MAX_MISMATCH = 2;
	
//STR parameters
int MAX_UNIT_LEN = 4;
int MIN_REPORT_UNITS = 3;
int MIN_REPORT_LEN = 7;
int DIST_FROM_STR = 1;

bool SCAFFOLD_CONTIGS = 0;
int INSERT_SIZE = 150;
int INSERT_STDEV = 15;

int num_windows = 0;

// constants
//////////////////////////////////////////////////////////////////////////

const char Graph_t::CODE_MAPPED = 'M';
const char Graph_t::CODE_BASTARD = 'B';

const string Graph_t::COLOR_ALL    = "white"; 
const string Graph_t::COLOR_LOW    = "grey";
const string Graph_t::COLOR_NOVO   = "darkorange3";
const string Graph_t::COLOR_TUMOR  = "red";
const string Graph_t::COLOR_NORMAL = "green";
const string Graph_t::COLOR_SHARED = "blue"; //"deepskyblue4";
const string Graph_t::COLOR_SOURCE = "orange\" style=\"filled";
const string Graph_t::COLOR_SINK   = "yellow\" style=\"filled";
const string Graph_t::COLOR_TOUCH  = "magenta";

/***********************************/

// print usage info to stderr
void printUsage();

// print help text to stderr
void printHelpText(Filters & filters);

// print configuration to file
void printConfiguration(ostream & out, Filters & filters);

// load referecne for fasta file
int loadRefs(const string reference, const string region, vector< map<string, Ref_t *> > &reftable, RefVector &bamrefs, int num_threads, int thread);

// loadbed : load regions from BED file
void loadBed(const string bedfile, vector< map<string, Ref_t *> > &reftable, RefVector &bamrefs, int num_threads);

static void* execute(void* ptr);

int rLancet(string tumor_bam, string normal_bam, string ref_fasta, string reg, string bed_file, int numthreads);

#endif
