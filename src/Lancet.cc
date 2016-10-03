#include "Microassembler.hh"

/****************************************************************************
** Lancet.cc
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

string VERSION = "1.0.0 (beta), September 18 2016";

/****  configuration parameters ****/
int NUM_THREADS = 1;

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
int WINDOW_SIZE = 600;

string TUMOR;
string NORMAL;
string RG_FILE;
string REFFILE;
string BEDFILE;
string REGION;

int minK = 11;
int maxK = 100;
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

void printConfiguration(ostream & out, Filters & filters)
{
	out << "tumor-BAM: "        << TUMOR << endl;
	out << "normal-BAM: "       << NORMAL << endl;
	out << "reference: "        << REFFILE << endl;
	out << "region: "           << REGION  << endl;
	out << "BED-file: "         << BEDFILE  << endl;

	out << "min-K: "            << minK << endl;
	out << "max-K: "            << maxK << endl;
	out << "tip-len: "          << MAX_TIP_LEN << endl;

	//out << "MIN_THREAD_READS: " << MIN_THREAD_READS << endl;
	out << "cov-thr: "          << COV_THRESHOLD << endl;
	cerr.unsetf(ios::floatfield); // floatfield not set
	cerr.precision(5);
	out << "cov-ratio: "        << MIN_COV_RATIO << endl;
	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	out << "low-cov: "          << LOW_COV_THRESHOLD << endl;
	out << "window-size: "      << WINDOW_SIZE << endl;
	out << "max-avg-cov: "      << MAX_AVG_COV << endl;
	out << "min-map-qual: "     << MIN_MAP_QUAL << endl;
	out << "min-base-qual: "    << MIN_QV_CALL << endl;
	out << "trim-lowqual: "     << MIN_QV_TRIM << endl;
	out << "quality-range: "    << QV_RANGE << endl;	
	out << "node-str-len: "     << NODE_STRLEN << endl;
	out << "dfs-limit: "        << DFS_LIMIT << endl;
	out << "max-indel-len: "    << MAX_INDEL_LEN << endl;
	out << "max-mismatch: "     << MAX_MISMATCH << endl;
	out << "num-threads: "      << NUM_THREADS << endl;	
	//out << "SCAFFOLD_CONTIGS: " << bvalue(SCAFFOLD_CONTIGS) << endl;
	//out << "INSERT_SIZE: "      << INSERT_SIZE << " +/- " << INSERT_STDEV << endl;
	
	//filters
	out << "min-phred-fisher: "     << filters.minPhredFisher << endl;
	out << "min-strand-bias: "      << filters.minStrandBias << endl;
	out << "min-alt-count-tumor: "  << filters.minAltCntTumor << endl;
	out << "max-alt-count-normal: " << filters.maxAltCntNormal << endl;
	out << "min-vaf-tumor: "        << filters.minVafTumor << endl;
	out << "max-vaf-normal: "       << filters.maxVafNormal << endl;
	out << "min-coverage-tumor: "   << filters.minCovTumor << endl;
	out << "max-coverage-tumor: "   << filters.maxCovTumor << endl;
	out << "min-coverage-normal: "  << filters.minCovNormal << endl;
	out << "max-coverage-normal: "  << filters.maxCovNormal << endl;
	
	out << "active-regions: "   << bvalue(ACTIVE_REGIONS) << endl;
	out << "kmer-recovery: "    << bvalue(KMER_RECOVERY) << endl;
	out << "print-graphs: "     << bvalue(PRINT_ALL) << endl;
	out << "verbose: "          << bvalue(verbose) << endl;
	out << "more-verbose: "     << bvalue(VERBOSE) << endl;
	
	out << endl;
}


// loadRef
//////////////////////////////////////////////////////////////
int loadRefs(const string reference, const string region, vector< map<string, Ref_t *> > &reftable, int num_threads, int thread)	
{
	//if(verbose) { cerr << "LoadRef " << reference << endl; }

	// open fasta index
    faidx_t *fai = fai_load(reference.c_str());
    if ( !fai ) { cerr << "Could not load fai index of " << reference << endl; }

	// extrat sequence
    int seq_len;
    char *seq = fai_fetch(fai, region.c_str(), &seq_len);
    if ( seq_len < 0 ) { cerr << "Failed to fetch sequence in " << region << endl; }
		
	// convert char* to string
	string s(seq, seq + seq_len);
	
	string ss;
	string hdr = region;

	// extrat coordinates for header
	size_t x     = hdr.find_first_of(':');
	size_t y     = hdr.find_first_of('-', x);
	string CHR   = hdr.substr(0,x);
	string START = hdr.substr(x+1, y-x-1);
	string END   = hdr.substr(y+1, string::npos);
	
	for (unsigned int i = 0; i < s.length(); ++i) { s[i] = toupper(s[i]); }

	// split into overalpping windows if sequence is too long
	int end = s.length();
	int offset = 0;
	int delta = 100;
	
	int T = thread; // thread counter
	for (; offset < end; offset+=delta) {
		
		// adjust end if 
		int LEN = WINDOW_SIZE;
		if( (offset + WINDOW_SIZE) > (int)s.length() ) { 
			LEN = s.length() - offset;
			end = offset;
		}
		
		ss = s.substr(offset,LEN);
		
		// make new reference entry 
		Ref_t * ref = new Ref_t(minK);
	
		ref->refchr   = CHR;
		ref->refstart = atoi(START.c_str()) + offset;
		ref->refend   = ref->refstart + LEN - 1;
		//ref->refend   = atoi(end.c_str());
		
		hdr = ref->refchr;
		hdr += ":";
		hdr += itos(ref->refstart);
		hdr += "-";
		hdr += itos(ref->refend);

		if(verbose) { cerr << "hdr:\t" << hdr << endl; }
		
		ref->setHdr(hdr);
		ref->setSeq(ss);
		ref->setRawSeq(ss);
		
		ref->hdr = hdr;
		
		reftable[T].insert(make_pair(hdr, ref));
		++num_windows;
		
		// move to next reftable
		++T;
		if( (num_windows%num_threads) == 0) { T=0; }
	}
	
	return T;
}

// loadbed : load regions from BED file
//////////////////////////////////////////////////////////////
void loadBed(const string bedfile, vector< map<string, Ref_t *> > &reftable, int num_threads) { 
	
	int num_regions = 0;
	string line;
	string region;
	vector<std::string> tokens;
	ifstream bfile (bedfile);
	int t = 0;
	if (bfile.is_open()) {
		while ( getline (bfile,line) ) {
			
			size_t x = line.find_first_of('#');
			if(x == 0) { continue; } // skip comments
			
			++num_regions;
			//cerr << line << '\n';
			
			// extrat coordinates			
		    istringstream iss(line);
		    string token;
			tokens.clear();
		    while(std::getline(iss, token, '\t')) {  // but we can specify a different one
				tokens.push_back(token);
			}

			region = tokens[0] + ":" + tokens[1] + "-" + tokens[2];			
			t = loadRefs(REFFILE,region,reftable,num_threads, t);
		}
		bfile.close();
		
		cerr << "Loaded " << num_regions << " from bedfile" << endl;
	}
	else {
		cerr << "Couldn't open " << bedfile << endl;
		exit(1);
	}
}

static void* execute(void* ptr) {

    Microassembler* ma = (Microassembler*)ptr;
	
	ma->processReads();
	//ma->vDB.printToVCF();
	
	pthread_exit(NULL);
}

// main
//////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{			
	stringstream LOGO1;
	LOGO1 <<
		"\n"
		"  |                           |  \n"
		"  |      _` | __ \\   __|  _ \\ __|\n"
		"  |     (   | |   | (     __/ |  \n"
		" _____|\\__,_|_|  _|\\___|\\___|\\__|\n";
	
	stringstream HEADER;
	HEADER << LOGO1.str() <<
		"\nProgram: lancet (micro-assembly somatic variant caller)\n"
		"Version: "<< VERSION << "\n"
		"Contact: Giuseppe Narzisi <gnarzisi@nygenome.org>\n";
	
	string USAGE = "\nUsage: lancet [options] --tumor <BAM file> --normal <BAM file> --ref <FASTA file> --reg <chr:start-end>\n [-h for full list of commands]\n\n";

	if (argc == 1)
	{
		cerr << HEADER.str() << USAGE;
		exit(0);
	}

	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	
	// initilize filter thresholds
	Filters filters;
	filters.minPhredFisher = 5;
	filters.minCovNormal = 10;
	filters.maxCovNormal = 1000000;
	filters.minCovTumor = 4;
	filters.maxCovTumor = 1000000;
	filters.minVafTumor = 0.05;
	filters.maxVafNormal = 0;
	filters.minAltCntTumor = 3;
	filters.maxAltCntNormal = 0;
	filters.minStrandBias = 1;
	
	stringstream helptext;
	helptext << HEADER.str() << USAGE <<
		"Required\n"
		"   --tumor, -t              <BAM file>    : BAM file of mapped reads for tumor\n"
		"   --normal, -n             <BAM file>    : BAM file of mapped reads for normal\n"
		"   --ref, -r                <FASTA file>  : FASTA file of reference genome\n"
		"   --reg, -p                <string>      : genomic region (in chr:start-end format)\n"
		"   --bed, -B                <string>      : genomic regions from file (BED format)\n"

		"\nOptional\n"
		"   --min-k, k                <int>         : min kmersize [default: " << minK << "]\n"
		"   --max-k, -K               <int>         : max kmersize [default: " << maxK << "]\n"
		"   --trim-lowqual, -q        <int>         : trim bases below qv at 5' and 3' [default: " << MIN_QV_TRIM << "]\n"
		"   --min-base-qual, -C       <int>         : minimum base quality required to consider a base for SNV calling [default: " << MIN_QV_CALL << "]\n"
		"   --quality-range, -Q       <char>        : quality value range [default: " << (char) QV_RANGE << "]\n"
		"   --min-map-qual, -b        <inr>         : minimum read mapping quality in Phred-scale [default: " << MIN_MAP_QUAL << "]\n"
		"   --tip-len, -l             <int>         : max tip length [default: " << MAX_TIP_LEN << "]\n"
		"   --cov-thr, -c             <int>         : coverage threshold [default: " << COV_THRESHOLD << "]\n"
		"   --cov-ratio, -x           <float>       : minimum coverage ratio [default: " << MIN_COV_RATIO << "]\n"
		"   --max-avg-cov, -u         <int>         : maximum average coverage allowed per region [default: " << MAX_AVG_COV << "]\n"
		"   --low-cov, -d             <int>         : low coverage threshold [default: " << LOW_COV_THRESHOLD << "]\n"
		"   --window-size, -w         <int>         : window size of the region to assemble (in base-pairs) [default: " << WINDOW_SIZE << "]\n"
		"   --dfs-limit, -F           <int>         : limit dfs/bfs graph traversal search space [default: " << DFS_LIMIT << "]\n"
		"   --max-indel-len, -T       <int>         : limit on size of detectable indel [default: " << MAX_INDEL_LEN << "]\n"
		"   --max-mismatch, -M        <int>         : max number of mismatches for near-perfect repeats [default: " << MAX_MISMATCH << "]\n"
		"   --num-threads, -X         <int>         : number of parallel threads [default: " << NUM_THREADS << "]\n"
//		"   --rg-file, -g             <string>      : read group file\n"
		"   --node-str-len, -L        <int>         : length of sequence to display at graph node (default: " << NODE_STRLEN << ")\n"

		"\nFilters\n"
		"   --min-alt-count-tumor, -a  <int>        : minimum alternative count in the tumor [default: " << filters.minAltCntTumor << "]\n"
		"   --max-alt-count-normal, -m <int>        : maximum alternative count in the normal [default: " << filters.maxAltCntNormal << "]\n"
		"   --min-vaf-tumor, -e        <float>      : minimum variant allele frequency (AlleleCov/TotCov) in the tumor [default: " << filters.minVafTumor << "]\n"
		"   --max-vaf-normal, -i       <float>      : maximum variant allele frequency (AlleleCov/TotCov) in the normal [default: " << filters.maxVafNormal << "]\n"
		"   --min-coverage-tumor, -o   <int>        : minimum coverage in the tumor [default: " << filters.minCovTumor << "]\n"
		"   --max-coverage-tumor, -y   <int>        : maximum coverage in the tumor [default: " << filters.maxCovTumor << "]\n"
		"   --min-coverage-normal, -z  <int>        : minimum coverage in the normal [default: " << filters.minCovNormal << "]\n"
		"   --max-coverage-normal, -j  <int>        : maximum coverage in the normal [default: " << filters.maxCovNormal << "]\n"
		"   --min-phred-fisher, -s     <float>      : minimum fisher exact test score [default: " << filters.minPhredFisher << "]\n"
		"   --min-strand-bias, -f      <float>      : minimum strand bias threshold [default: " << filters.minStrandBias << "]\n"
			
		"\nShort Tandem Repeat parameters\n"
		"   --max-unit-length, -U      <int>        : maximum unit length of the motif [default: " << MAX_UNIT_LEN << "]\n"
		"   --min-report-unit, -N      <int>        : minimum number of units to report [default: " << MIN_REPORT_UNITS << "]\n"
		"   --min-report-len, -Y       <int>        : minimum length of tandem in base pairs [default: " << MIN_REPORT_LEN << "]\n"
		"   --dist-from-str, -D        <int>        : distance (in bp) of variant from STR locus [default: " << DIST_FROM_STR << "]\n"
		
		"\nFlags\n"
		"   --active-region-off, -W    : turn off active region module\n"		
		"   --kmer-recovery, -R        : turn on k-mer recovery (experimental)\n"
		"   --print-graph, -A          : print graph (in .dot format) after every stage\n"
		"   --verbose, -v              : be verbose\n"
		"   --more-verbose, -V         : be more verbose\n"
		"\n";

	bool errflg = false;
	int ch;

	optarg = NULL;
	
	static struct option long_options[] = {
		
		// required
		{"tumor",   required_argument, 0, 't'},
		{"normal",  required_argument, 0, 'n'},
		{"ref",     required_argument, 0, 'r'},
		{"bed",     required_argument, 0, 'B'},
		
		// optional
		{"reg",  required_argument, 0, 'p'},
		{"rg-file",  required_argument, 0, 'g'},
		{"min-k",    required_argument, 0, 'k'},
		{"max-k",    required_argument, 0, 'K'},
		{"tip-len",  required_argument, 0, 'l'},
		{"cov-thr",  required_argument, 0, 'c'},
		{"cov-ratio",  required_argument, 0, 'x'},
		{"low-cov",  required_argument, 0, 'd'},
		{"window-size",  required_argument, 0, 'w'},
		{"max-avg-cov",  required_argument, 0, 'u'},
		{"min-map-qual",  required_argument, 0, 'b'},
		{"min-base-qual",  required_argument, 0, 'C'},
		{"trim-lowqual",  required_argument, 0, 'q'},
		{"quality-range",  required_argument, 0, 'Q'},

		{"node-str-len",  required_argument, 0, 'L'},
		{"dfs-limit",  required_argument, 0, 'F'},
		//{"path-limit",  required_argument, 0, 'P'},
		{"num-threads",  required_argument, 0, 'X'},
		{"max-indel-len",  required_argument, 0, 'T'},
		{"max-mismatch",  required_argument, 0, 'M'},

		// STR params
		{"max-unit-length",  required_argument, 0, 'U'},
		{"min-report-unit",  required_argument, 0, 'N'},
		{"min-report-len",  required_argument, 0, 'Y'},
		{"dist-from-str",  required_argument, 0, 'D'},
		
		//filters
		{"min-phred-fisher",  required_argument, 0, 's'},
		{"min-strand-bias",  required_argument, 0, 'f'},
		{"min-alt-count-tumor",  required_argument, 0, 'a'},
		{"max-alt-count-normal",  required_argument, 0, 'm'},
		{"min-vaf-tumor",  required_argument, 0, 'e'},
		{"max-vaf-normal",  required_argument, 0, 'i'},
		{"min-coverage-tumor",  required_argument, 0, 'o'},
		{"max-coverage-tumor",  required_argument, 0, 'y'},
		{"min-coverage-normal",  required_argument, 0, 'z'},
		{"max-coverage-normal",  required_argument, 0, 'j'},

		{"active-region-off", no_argument,      0, 'W'},		
		{"kmer-recovery-on", no_argument,      0, 'R'},		
		{"erroflag", no_argument,      0, 'h'},		
		{"verbose", no_argument,       0, 'v'},
		{"more-verbose", no_argument,  0, 'V'},
		{"print-graph", no_argument,     0, 'A'},
		
		{0, 0, 0, 0}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;

	//while (!errflg && ((ch = getopt (argc, argv, "u:m:n:r:g:s:k:K:l:t:c:d:x:BDRACIhSL:T:M:vF:q:b:Q:P:p:E")) != EOF))
	while (!errflg && ((ch = getopt_long (argc, argv, "u:n:r:g:k:K:l:f:t:c:C:d:x:ARhSWL:T:M:vVF:q:b:B:Q:p:s:a:m:e:i:o:y:z:w:j:X:U:N:Y:D:", long_options, &option_index)) != -1))
	{
		switch (ch)
		{
			case 't': TUMOR            = optarg;       break; 
			case 'n': NORMAL           = optarg;       break; 
			case 'r': REFFILE          = optarg;       break;
			case 'B': BEDFILE          = optarg;       break;
			case 'p': REGION           = optarg;       break;
			
			case 'g': RG_FILE          = optarg;       break;
			case 'k': minK             = atoi(optarg); break;
			case 'K': maxK             = atoi(optarg); break;
			case 'l': MAX_TIP_LEN      = atoi(optarg); break;
			//case 't': MIN_THREAD_READS = atoi(optarg); break;
			case 'c': COV_THRESHOLD    = atoi(optarg); break;
			case 'x': MIN_COV_RATIO    = atof(optarg); break;
			case 'd': LOW_COV_THRESHOLD= atoi(optarg); break;
			case 'w': WINDOW_SIZE      = atoi(optarg); break;
			case 'u': MAX_AVG_COV      = atoi(optarg); break;
			
			case 'q': MIN_QV_TRIM      = atoi(optarg); break;
			case 'C': MIN_QV_CALL      = atoi(optarg); break;
			case 'b': MIN_MAP_QUAL     = atoi(optarg); break;
			case 'Q': QV_RANGE         = *optarg;      break;

			case 'L': NODE_STRLEN      = atoi(optarg); break;
			case 'F': DFS_LIMIT        = atoi(optarg); break;
			case 'X': NUM_THREADS      = atoi(optarg); break;
			case 'T': MAX_INDEL_LEN    = atoi(optarg); break;
			case 'M': MAX_MISMATCH     = atoi(optarg); break;
			
			case 's': filters.minPhredFisher = atoi(optarg); break;
			case 'f': filters.minStrandBias = atoi(optarg); break;
			case 'a': filters.minAltCntTumor = atoi(optarg); break;
			case 'm': filters.maxAltCntNormal = atoi(optarg); break;
			case 'e': filters.minVafTumor = atoi(optarg); break;
			case 'i': filters.maxVafNormal = atoi(optarg); break;
			case 'o': filters.minCovTumor = atoi(optarg); break;
			case 'y': filters.maxCovTumor = atoi(optarg); break;
			case 'z': filters.minCovNormal = atoi(optarg); break;
			case 'j': filters.maxCovNormal = atoi(optarg); break;

			case 'W': ACTIVE_REGIONS   = 0;            break;
			case 'R': KMER_RECOVERY    = 1;            break;
			case 'v': verbose          = 1;            break;
			case 'V': VERBOSE=1; verbose=1;            break;
			case 'A': PRINT_ALL        = 1;            break;

			case 'h': errflg = 1;                      break;

			case '?':
			fprintf (stderr, "Unrecognized option -%c\n", optopt);

			default:
			errflg = true;
		}

		if (errflg)
		{
			cout << helptext.str();
			exit (EXIT_FAILURE);
		}
	}

	// update min base quality values
	MIN_QUAL_TRIM = MIN_QV_TRIM + QV_RANGE;
	MIN_QUAL_CALL = MIN_QV_CALL + QV_RANGE;

	if (TUMOR == "") { cerr << "ERROR: Must provide the tumor BAM file (-t)" << endl; ++errflg; }
	if (NORMAL == "") { cerr << "ERROR: Must provide the normal BAM file (-n)" << endl; ++errflg; }		
	if (REFFILE == "") { cerr << "ERROR: Must provide a reference genome file (-r)" << endl; ++errflg; }
	if ( (BEDFILE == "") && (REGION == "") ) { cerr << "ERROR: Must provide region (-p) or BED file (-B)" << endl; ++errflg; }

	if (errflg) { exit(EXIT_FAILURE); }
	
    ofstream params_file;
    params_file.open ("config.txt");
	printConfiguration(params_file, filters); // save parameters setting to file
    params_file.close();

	if(verbose) { printConfiguration(cerr, filters); }
	
	// run the assembler on each region
	try {
		
		//Microassembler* driver = new Microassembler();
		//driver->run(argc, argv);
		//assembler->vDB.printToVCF();

		pthread_t threads[NUM_THREADS];
		pthread_attr_t attr;
		void * status;
		int rc;
		int i;		
		vector<Microassembler*> assemblers(NUM_THREADS, new Microassembler());
		vector< map<string, Ref_t *> > reftables(NUM_THREADS, map<string, Ref_t *>()); // table of references to analyze
		
		if (BEDFILE != "") {
			loadBed(BEDFILE,reftables,NUM_THREADS);
		}
		if (REGION != "") {
			loadRefs(REFFILE,REGION,reftables,NUM_THREADS, 0);
		}
		
		cerr << num_windows << " total windows to process" << endl << endl;
		
		// Initialize and set thread joinable
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		for( i=0; i < NUM_THREADS; ++i ) {
			cerr << "starting thread " << (i+1) << " on " << reftables[i].size() << " windows" << endl;
		
			assemblers[i] = new Microassembler();

			assemblers[i]->ACTIVE_REGION_MODULE = ACTIVE_REGIONS;
			assemblers[i]->KMER_RECOVERY = KMER_RECOVERY;
			assemblers[i]->verbose = verbose;
			assemblers[i]->VERBOSE = VERBOSE;
			assemblers[i]->PRINT_DOT_READS = PRINT_DOT_READS;
			assemblers[i]->PRINT_ALL = PRINT_ALL;
			assemblers[i]->MIN_QV_CALL = MIN_QV_CALL;
			assemblers[i]->MIN_QV_TRIM = MIN_QV_TRIM;
			assemblers[i]->QV_RANGE = QV_RANGE;
			assemblers[i]->MIN_QUAL_TRIM = MIN_QUAL_TRIM;
			assemblers[i]->MIN_QUAL_CALL = MIN_QUAL_CALL;
			assemblers[i]->MIN_MAP_QUAL = MIN_MAP_QUAL;
			assemblers[i]->TUMOR = TUMOR;
			assemblers[i]->NORMAL = NORMAL;
			assemblers[i]->RG_FILE = RG_FILE;
			assemblers[i]->REFFILE = REFFILE;
			assemblers[i]->minK = minK;
			assemblers[i]->maxK = maxK;
			assemblers[i]->MAX_TIP_LEN = MAX_TIP_LEN;
			assemblers[i]->MIN_THREAD_READS = MIN_THREAD_READS;
			assemblers[i]->COV_THRESHOLD = COV_THRESHOLD;
			assemblers[i]->MIN_COV_RATIO = MIN_COV_RATIO;
			assemblers[i]->LOW_COV_THRESHOLD = LOW_COV_THRESHOLD;
			assemblers[i]->MAX_AVG_COV = MAX_AVG_COV;
			assemblers[i]->NODE_STRLEN = NODE_STRLEN;
			assemblers[i]->DFS_LIMIT = DFS_LIMIT;
			assemblers[i]->MAX_INDEL_LEN = MAX_INDEL_LEN;
			assemblers[i]->MAX_MISMATCH = MAX_MISMATCH;		
			assemblers[i]->MAX_UNIT_LEN = MAX_UNIT_LEN;
			assemblers[i]->MIN_REPORT_UNITS = MIN_REPORT_UNITS;
			assemblers[i]->MIN_REPORT_LEN = MIN_REPORT_LEN;
			assemblers[i]->DIST_FROM_STR = DIST_FROM_STR;	
			
			assemblers[i]->reftable = &reftables[i];
			assemblers[i]->setFilters(&filters);
			assemblers[i]->setID(i+1);
	
			rc = pthread_create(&threads[i], NULL, execute, (void * )assemblers[i]);
			
			if (rc){
				cerr << "Error:unable to create thread," << rc << endl;
				exit(-1);
			}
		}
		
		// free attribute and wait for the other threads
		pthread_attr_destroy(&attr);
		for( i=0; i < NUM_THREADS; ++i ){
			rc = pthread_join(threads[i], &status);
			if (rc){
				cerr << "Error:unable to join," << rc << endl;
				exit(-1);
			}
			cerr << "Main: completed thread id :" << (i+1) ;
			cerr << " exiting with status :" << status << endl;
		}
		
		int tot_skip = 0;
		int tot_svn_only = 0;
		int tot_indel_only = 0;
		int tot_softclip_only = 0;
		int tot_indel_or_softclip = 0;
		int tot_snv_or_indel = 0;
		int tot_snv_or_softclip = 0;
		int tot_snv_or_indel_or_softclip = 0;
		//merge variant from all threads
		cerr << "Merge variants" << endl;
		VariantDB_t variantDB; // variants DB
		for( i=0; i < NUM_THREADS; ++i ) {
			
			tot_skip += assemblers[i]->num_skip;
			tot_svn_only += assemblers[i]->num_snv_only_regions;
			tot_indel_only += assemblers[i]->num_indel_only_regions;
			tot_softclip_only += assemblers[i]->num_softclip_only_regions;
			tot_indel_or_softclip += assemblers[i]->num_indel_or_softclip_regions;			
			tot_snv_or_indel += assemblers[i]->num_snv_or_indel_regions;			
			tot_snv_or_softclip += assemblers[i]->num_snv_or_softclip_regions;
			tot_snv_or_indel_or_softclip += assemblers[i]->num_snv_or_indel_or_softclip_regions;
						
			map<string,Variant_t> db = (assemblers[i]->vDB).DB;
			map<string,Variant_t>::iterator it;			
			for (it=db.begin(); it!=db.end(); ++it) {
				variantDB.addVar(it->second);
			}
		}
		
		//if(verbose) {
			cerr << "Total # of skipped windows: " << tot_skip << " (" << (100*(double)tot_skip/double(num_windows)) << "\%)" << endl;
			cerr << "- # of windows with SNVs only: " << tot_svn_only << endl;
			cerr << "- # of windows with indels only: " << tot_indel_only << endl;
			cerr << "- # of windows with softclips only: " << tot_softclip_only << endl;
			cerr << "- # of windows with indels or softclips: " << tot_indel_or_softclip << endl;
			cerr << "- # of windows with SNVs or indels: " << tot_snv_or_indel << endl;
			cerr << "- # of windows with SNVs or softclips: " << tot_snv_or_softclip << endl;
			cerr << "- # of windows with SNVs or indels or softclips: " << tot_snv_or_indel_or_softclip << endl;
		//}
		
		/***** get current time and date *****/
		time_t rawtime;
		time (&rawtime);
		char* DATE = ctime (&rawtime);
		/***************************************/
		
		variantDB.printToVCF(VERSION, REFFILE, DATE, filters, assemblers[0]->sample_name_normal, assemblers[0]->sample_name_tumor);		
	}
	catch (int e) {
		cerr << "An exception occurred. Exception Nr. " << e << endl;
	}

	pthread_exit(NULL);
}
