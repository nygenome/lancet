#include "Microassembler.hh"

/******************************************************************
** Lancet.cc
**
** Tool for localized assembly of genomic regions using
** the de Bruijn paradigm to detect genetic variants
** 
**  Authors: Giuseppe Narzisi 
**    Date: January 16, 2016
**
*******************************************************************/

string VERSION = "1.0.0 (beta), January 27 2016";

/****  configuration parameters ****/
int NUM_THREADS = 1;

bool verbose = false;
bool VERBOSE = false;
bool PRINT_ALL = false;
bool PRINT_DOT_READS = true;
int MIN_QV = 10;
int QV_RANGE = '!';
int MIN_QUAL = MIN_QV + QV_RANGE;
int MIN_MAP_QUAL = 0;
int WINDOW_SIZE = 600;

string TUMOR;
string NORMAL;
string RG_FILE;
string REFFILE;
string REGION;

int minK = 11;
int maxK = 100;
int MAX_TIP_LEN = minK;
unsigned int MIN_THREAD_READS = 3;
int COV_THRESHOLD = 5;
double MIN_COV_RATIO = 0.01;
int LOW_COV_THRESHOLD = 1;
int MAX_AVG_COV = 10000;
int NODE_STRLEN = 100;
int DFS_LIMIT = 1000000;
int PATH_LIMIT = 0;
int MAX_INDEL_LEN = 250;
int MAX_MISMATCH = 3;

bool SCAFFOLD_CONTIGS = 0;
int INSERT_SIZE = 150;
int INSERT_STDEV = 15;

// constants
//////////////////////////////////////////////////////////////////////////

const char Graph_t::CODE_MAPPED = 'M';
const char Graph_t::CODE_BASTARD = 'B';

const string Graph_t::COLOR_ALL    = "white"; 
const string Graph_t::COLOR_LOW    = "grey";
const string Graph_t::COLOR_NOVO   = "darkorange3";
const string Graph_t::COLOR_TUMOR  = "red";
const string Graph_t::COLOR_NORMAL = "green";
const string Graph_t::COLOR_SHARED = "deepskyblue4";
const string Graph_t::COLOR_SOURCE = "orange\" style=\"filled";
const string Graph_t::COLOR_SINK   = "yellow\" style=\"filled";
const string Graph_t::COLOR_TOUCH  = "magenta";

/***********************************/

void printConfiguration(ostream & out)
{
	out << "VERBOSE: "          << verbose << endl;
	out << "MORE VERBOSE: "     << VERBOSE << endl;

	out << "tumor BAM: "        << TUMOR << endl;
	out << "normal BAM: "       << NORMAL << endl;
	out << "reffile: "          << REFFILE << endl;
	out << "region: "           << REGION  << endl;

	out << "min K: "            << minK << endl;
	out << "max K: "            << maxK << endl;
	out << "MAX_TIP_LEN: "      << MAX_TIP_LEN << endl;

	out << "QV_RANGE: "         << QV_RANGE << endl;
	out << "MIN_QV: "           << MIN_QV << endl;
	out << "MIN_QUAL: "         << (char) MIN_QUAL << endl;
	out << "MIN_MAP_QUAL: "     << MIN_MAP_QUAL << endl;
	out << "MAX_AVG_COV: "		<< MAX_AVG_COV << endl;

	out << "MIN_THREAD_READS: " << MIN_THREAD_READS << endl;
	out << "COV_THRESHOLD: "    << COV_THRESHOLD << endl;
	cerr.unsetf(ios::floatfield); // floatfield not set
	cerr.precision(5);
	out << "MIN_COV_RATIO: "    << MIN_COV_RATIO << endl;
	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	out << "LOW_COV_THRESHOLD: "<< LOW_COV_THRESHOLD << endl;
	out << "WINDOW_SIZE: "      << WINDOW_SIZE << endl;
	out << "DFS_LIMIT: "        << DFS_LIMIT << endl;
	out << "PATH_LIMIT: "       << PATH_LIMIT << endl;
	out << "MAX_INDEL_LEN: "    << MAX_INDEL_LEN << endl;
	out << "MAX_MISMATCH: "     << MAX_MISMATCH << endl;

	out << "SCAFFOLD_CONTIGS: " << bvalue(SCAFFOLD_CONTIGS) << endl;
	out << "INSERT_SIZE: "      << INSERT_SIZE << " +/- " << INSERT_STDEV << endl;

	out << "PRINT_ALL: "        << bvalue(PRINT_ALL) << endl;

	out << "NODE_STRLEN: "      << NODE_STRLEN << endl;
	
	out << endl;
}

// loadRef
//////////////////////////////////////////////////////////////

//void loadRefs(const string & filename, map<string, Ref_t *> &reftable, int NUM_THREADS)
void loadRefs(const string reference, const string region, vector< map<string, Ref_t *> > &reftable, int num_threads)	
{
	if(verbose) { cerr << "LoadRef " << reference << endl; }

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
	
	for (unsigned int i = 0; i < s.length(); i++) { s[i] = toupper(s[i]); }

	// split into overalpping windoes if sequence is too long
	int end = s.length();
	int offset = 0;
	int delta = 100;
	
	int T = 0; // thread counter
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
		
		// move to next reftable
		T++;
		if( (T%num_threads) == 0) { T=0; }
	}


	if(verbose) { cerr << "Loaded " << reftable.size() << " ref sequences" << endl << endl; }
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
	stringstream HEADER;
	HEADER << 
		"\nProgram: Lancet (micro-assembly somatic variant caller)\n"
		"Version: "<< VERSION << "\n"
		"Contact: Giuseppe Narzisi <gnarzisi@nygenome.org>\n";
	
	string USAGE = "\nUsage: Lancet [options] --tumor <BAM file> --normal <BAM file> --ref <FASTA file> --reg <chr:start-end>\n [-h for full list of commands]\n\n";

	if (argc == 1)
	{
		cerr << HEADER.str() << USAGE;
		exit(0);
	}

	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	
	// initilize filter thresholds
	Filters filters;
	filters.minPhredFisher = 10;
	filters.minCovNormal = 10;
	filters.maxCovNormal = 1000000;
	filters.minCovTumor = 4;
	filters.maxCovTumor = 1000000;
	filters.minVafTumor = 0.05;
	filters.maxVafNormal = 0;
	filters.minAltCntTumor = 4;
	filters.maxAltCntNormal = 0;
	
	stringstream helptext;
	helptext << HEADER.str() << USAGE <<
		"Required\n"
		"   --tumor, -t              <BAM file>    : BAM file of mapped reads for tumor\n"
		"   --normal, -n             <BAM file>    : BAM file of mapped reads for normal\n"
		"   --ref, -r                <FASTA file>  : FASTA file of reference genome\n"
		"   --reg, -p                <string>      : genomic region (in chr:start-end format)\n"
		"\nOptional\n"
		"   --min-k, k                <int>         : min kmersize [default: " << minK << "]\n"
		"   --max-k, -K               <int>         : max kmersize [default: " << maxK << "]\n"
		"   --trim-lowqual, -q        <int>         : trim bases below qv at 5' and 3' [default: " << MIN_QV << "]\n"
		"   --quality-range, -Q       <char>        : quality value range [default: " << (char) QV_RANGE << "]\n"
		"   --min-map-qual, -b        <inr>         : minimum read mapping quality in Phred-scale [default: " << MIN_MAP_QUAL << "]\n"
		"   --tip-len, -l             <int>         : max tip length [default: " << MAX_TIP_LEN << "]\n"
		"   --cov-thr, -c             <int>         : coverage threshold [default: " << COV_THRESHOLD << "]\n"
		"   --cov-ratio, -x           <float>       : minimum coverage ratio [default: " << MIN_COV_RATIO << "]\n"
		"   --max-avg-cov, -u         <int>         : maximum average coverage allowed per region [default: " << MAX_AVG_COV << "]\n"
		"   --low-cov, -d             <int>         : low coverage threshold [default: " << LOW_COV_THRESHOLD << "]\n"
		"   --window-size, -w         <int>         : window size of the region to assemble (in base-pairs) [default: " << WINDOW_SIZE << "]\n"
		"   --dfs-limit, -F           <int>         : limit dfs search space [default: " << DFS_LIMIT << "]\n"
		"   --path-limit, -P          <int>         : limit on number of paths to report [default: " << PATH_LIMIT << "]\n"
		"   --max-indel-len, -T       <int>         : limit on size of detectable indel [default: " << MAX_INDEL_LEN << "]\n"
		"   --max-mismatch, -M        <int>         : max number of mismatches for near-perfect repeats [default: " << MAX_MISMATCH << "]\n"
		"   --num-threads, -X         <int>         : number of parallel threads [default: " << NUM_THREADS << "]\n"
		"   --rg-file, -g             <string>      : read group file\n"

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
		
		"\nFlags\n"
		"   -A            : print graph (in .dot format) after every stage\n"
		"   -L <len>      : length of sequence to display at graph node (default: " << NODE_STRLEN << ")\n"
		"   -v            : be verbose\n"
		"   -V            : be more verbose\n"
		"\n";

	bool errflg = false;
	int ch;

	optarg = NULL;
	
	static struct option long_options[] = {
		
		// required
		{"tumor",   required_argument, 0, 't'},
		{"normal",  required_argument, 0, 'n'},
		{"ref",     required_argument, 0, 'r'},
		
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
		{"trim-lowqual",  required_argument, 0, 'q'},
		{"quality-range",  required_argument, 0, 'Q'},

		{"node-str-len",  required_argument, 0, 'L'},
		{"dfs-limit",  required_argument, 0, 'F'},
		{"path-limit",  required_argument, 0, 'P'},
		{"num-threads",  required_argument, 0, 'X'},
		{"max-indel-len",  required_argument, 0, 'T'},
		{"max-mismatch",  required_argument, 0, 'M'},
		
		//filters
		{"min-phred-fisher",  required_argument, 0, 's'},
		{"min-alt-count-tumor",  required_argument, 0, 'a'},
		{"max-alt-count-normal",  required_argument, 0, 'm'},
		{"min-vaf-tumor",  required_argument, 0, 'e'},
		{"max-vaf-normal",  required_argument, 0, 'i'},
		{"min-coverage-tumor",  required_argument, 0, 'o'},
		{"max-coverage-tumor",  required_argument, 0, 'y'},
		{"min-coverage-normal",  required_argument, 0, 'z'},
		{"max-coverage-normal",  required_argument, 0, 'j'},

		{"erroflag", no_argument,      0, 'h'},		
		{"verbose", no_argument,       0, 'v'},
		{"more-verbose", no_argument,  0, 'V'},
		{"print-denovo", no_argument,  0, 'D'},
		{"print-refpath", no_argument, 0, 'R'},
		{"print-all", no_argument,     0, 'A'},
		{"print-raw", no_argument,     0, 'I'},
		{"include-unmapped", no_argument, 0, 'B'},
		{"input-bam", no_argument,     0, 'C'},
		{"input-fastq", no_argument,   0, 'E'},
		
		{0, 0, 0, 0}
	};
	/* getopt_long stores the option index here. */
	int option_index = 0;

	//while (!errflg && ((ch = getopt (argc, argv, "u:m:n:r:g:s:k:K:l:t:c:d:x:BDRACIhSL:T:M:vF:q:b:Q:P:p:E")) != EOF))
	while (!errflg && ((ch = getopt_long (argc, argv, "u:n:r:g:k:K:l:t:c:d:x:AhSL:T:M:vVF:q:b:Q:P:p:s:a:m:e:i:o:y:z:w:j:X:", long_options, &option_index)) != -1))
	{
		switch (ch)
		{
			case 't': TUMOR            = optarg;       break; 
			case 'n': NORMAL           = optarg;       break; 
			case 'r': REFFILE          = optarg;       break;
			
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
			
			case 'q': MIN_QV           = atoi(optarg); break;
			case 'b': MIN_MAP_QUAL     = atoi(optarg); break;
			case 'Q': QV_RANGE         = *optarg;      break;

			case 'L': NODE_STRLEN      = atoi(optarg); break;
			case 'F': DFS_LIMIT        = atoi(optarg); break;
			case 'P': PATH_LIMIT       = atoi(optarg); break;
			case 'X': NUM_THREADS      = atoi(optarg); break;
			case 'T': MAX_INDEL_LEN    = atoi(optarg); break;
			case 'M': MAX_MISMATCH     = atoi(optarg); break;
			
			case 's': filters.minPhredFisher = atoi(optarg); break;
			case 'a': filters.minAltCntTumor = atoi(optarg); break;
			case 'm': filters.maxAltCntNormal = atoi(optarg); break;
			case 'e': filters.minVafTumor = atoi(optarg); break;
			case 'i': filters.maxVafNormal = atoi(optarg); break;
			case 'o': filters.minCovTumor = atoi(optarg); break;
			case 'y': filters.maxCovTumor = atoi(optarg); break;
			case 'z': filters.minCovNormal = atoi(optarg); break;
			case 'j': filters.maxCovNormal = atoi(optarg); break;

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

	if (TUMOR == "") { cerr << "ERROR: Must provide the tumor BAM file (-t)" << endl; errflg++; }
	if (NORMAL == "") { cerr << "ERROR: Must provide the normal BAM file (-n)" << endl; errflg++; }		
	if (REFFILE == "") { cerr << "ERROR: Must provide a reference genome file (-r)" << endl; errflg++; }

	if (errflg) { exit(EXIT_FAILURE); }

	if(verbose) { printConfiguration(cerr); }
	
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
		
		loadRefs(REFFILE,REGION,reftables,NUM_THREADS);

		// Initialize and set thread joinable
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		for( i=0; i < NUM_THREADS; i++ ) {
			cerr << "starting thread " << (i+1) << " on " << reftables[i].size() << " windows" << endl;
		
			assemblers[i] = new Microassembler();
						
			assemblers[i]->verbose = verbose;
			assemblers[i]->VERBOSE = VERBOSE;
			assemblers[i]->PRINT_DOT_READS = PRINT_DOT_READS;
			assemblers[i]->PRINT_ALL = PRINT_ALL;
			assemblers[i]->MIN_QV = MIN_QV;
			assemblers[i]->QV_RANGE = QV_RANGE;
			assemblers[i]->MIN_QUAL = MIN_QUAL;
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
			assemblers[i]->PATH_LIMIT = PATH_LIMIT;
			assemblers[i]->MAX_INDEL_LEN = MAX_INDEL_LEN;
			assemblers[i]->MAX_MISMATCH = MAX_MISMATCH;	
			
			assemblers[i]->reftable = &reftables[i];
			assemblers[i]->setFilters(filters);
			assemblers[i]->setID(i+1);
	
			rc = pthread_create(&threads[i], NULL, execute, (void * )assemblers[i]);
			
			if (rc){
				cerr << "Error:unable to create thread," << rc << endl;
				exit(-1);
			}
		}
		
		// free attribute and wait for the other threads
		pthread_attr_destroy(&attr);
		for( i=0; i < NUM_THREADS; i++ ){
			rc = pthread_join(threads[i], &status);
			if (rc){
				cerr << "Error:unable to join," << rc << endl;
				exit(-1);
			}
			cerr << "Main: completed thread id :" << (i+1) ;
			cerr << " exiting with status :" << status << endl;
		}
		
		//merge variant from all threads
		cerr << "Merge variants" << endl;
		VariantDB_t variantDB; // variants DB
		for( i=0; i < NUM_THREADS; i++ ) {		
			map<string,Variant_t> db = (assemblers[i]->vDB).DB;
			map<string,Variant_t>::iterator it;			
			for (it=db.begin(); it!=db.end(); ++it) {
				variantDB.addVar(it->second);
			}
		}
		
		
		/***** get current time and date *****/
		
		time_t rawtime;
		time (&rawtime);
		char* DATE = ctime (&rawtime);
		
		/***************************************/
		
		variantDB.printToVCF(VERSION,REFFILE,DATE,filters, assemblers[0]->sample_name_normal, assemblers[0]->sample_name_tumor);
	}
	catch (int e) {
		cerr << "An exception occurred. Exception Nr. " << e << endl;
	}

	pthread_exit(NULL);
}
