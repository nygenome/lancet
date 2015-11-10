#include "Lancet.hh"

/******************************************************************
** Microassembler.cc
**
** Tool for localized assembly of genomic regions using
** the de Bruijn paradigm to detect genetic variants
** 
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

void Microassembler::printConfiguration(ostream & out)
{
	out << "VERBOSE: "          << VERBOSE << endl;

	out << "tumor BAM: "        << TUMOR << endl;
	out << "normal BAM: "       << NORMAL << endl;
	out << "reffile: "          << REFFILE << endl;
	out << "prefix: "           << PREFIX  << endl;

	out << "min K: "            << minK << endl;
	out << "max K: "            << maxK << endl;
	out << "MAX_TIP_LEN: "      << MAX_TIP_LEN << endl;

	out << "QV_RANGE: "         << QV_RANGE << endl;
	out << "MIN_QV: "           << MIN_QV << endl;
	out << "MIN_QUAL: "         << (char) MIN_QUAL << endl;
	out << "MIN_MAP_QUAL: "     << MIN_MAP_QUAL << endl;
	out << "MAX_AVG_COV: "		<< MAX_AVG_COV << endl;

	out << "INCLUDE_BASTARDS: " << bvalue(INCLUDE_BASTARDS) << endl;

	out << "MIN_THREAD_READS: " << MIN_THREAD_READS << endl;
	out << "COV_THRESHOLD: "    << COV_THRESHOLD << endl;
	cerr.unsetf(ios::floatfield); // floatfield not set
	cerr.precision(5);
	out << "MIN_COV_RATIO: "    << MIN_COV_RATIO << endl;
	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);
	out << "TIP_COV_THRESHOLD: "<< TIP_COV_THRESHOLD << endl;
	out << "DFS_LIMIT: "        << DFS_LIMIT << endl;
	out << "PATH_LIMIT: "       << PATH_LIMIT << endl;
	out << "MAX_INDEL_LEN: "    << MAX_INDEL_LEN << endl;
	out << "MAX_MISMATCH: "     << MAX_MISMATCH << endl;

	out << "SCAFFOLD_CONTIGS: " << bvalue(SCAFFOLD_CONTIGS) << endl;
	out << "INSERT_SIZE: "      << INSERT_SIZE << " +/- " << INSERT_STDEV << endl;

	out << "PRINT_ALL: "        << bvalue(PRINT_ALL) << endl;
	out << "PRINT_RAW: "        << bvalue(PRINT_RAW) << endl;

	out << "PRINT_DENOVO: "     << bvalue(PRINT_DENOVO) << endl;
	out << "PRINT_REFPATH: "    << bvalue(PRINT_REFPATH) << endl;
	out << "NODE_STRLEN: "      << NODE_STRLEN << endl;

    out << "FASTQ_ASM: "        << bvalue(FASTQ_ASM) << endl;
}

// loadRef
//////////////////////////////////////////////////////////////

void Microassembler::loadRefs(const string & filename)
{
	cerr << "LoadRef " << filename << endl;

	FILE * fp = xfopen(filename, "r");

	string s, hdr;

	while (Fasta_Read (fp, s, hdr))
	{
		for (unsigned int i = 0; i < s.length(); i++) { s[i] = toupper(s[i]); }

		Ref_t * ref = new Ref_t(minK);
		
		ref->setHdr(hdr);
		ref->setSeq(s);
		ref->setRawSeq(s);

		cerr << "hdr:\t" << hdr << endl;

		size_t x = hdr.find_first_of(':');
		size_t y = hdr.find_first_of('-', x);

		string start = hdr.substr(x+1, y-x-1);
		string end   = hdr.substr(y+1, string::npos);

		ref->refchr   = hdr.substr(0,x);
		ref->refstart = atoi(start.c_str());
		ref->refend   = atoi(end.c_str());

		size_t z = hdr.find_first_of(';', y);

		if (z != string::npos)
		{
			y = hdr.find_first_of('-', z);

			start = hdr.substr(z+1, y-z-1);
			end =   hdr.substr(y+1, string::npos);

			hdr = ref->refchr;
			hdr += ":";
			hdr += start;
			hdr += "-";
			hdr += end;

			ref->hdr = hdr;
		}

		//cerr << "label:\t"    << ref->hdr << endl;
		//cerr << "refchr:\t"   << ref->refchr << endl;
		//cerr << "refstart:\t" << ref->refstart << endl;
		//cerr << "refend:\t"   << ref->refend << endl;

		reftable.insert(make_pair(hdr, ref));
	}

	cerr << "Loaded " << reftable.size() << " ref sequences" << endl << endl;

	xfclose(fp);
}

// load Red Groups
//////////////////////////////////////////////////////////////

void Microassembler::loadRG(const string & filename, int member) {

	FILE * fp;
	
	switch (member) {
		case FATHER:
			fp = xfopen(PREFIX + "/father.rg", "r");
			break;
		case MOTHER:
			fp = xfopen(PREFIX + "/mother.rg", "r");
			break;
		case SELF:
			fp = xfopen(PREFIX + "/self.rg", "r");
			break;
		case SIBLING:
			fp = xfopen(PREFIX + "/sibling.rg", "r");
			break;
		default:
			fp = xfopen(PREFIX + "/" + filename, "r");
	}
	
	char rgbuffer[BUFFER_SIZE];
	string rg;

	while (fscanf(fp, "%s\n", rgbuffer) == 1) {
		//memcpy(a,s.c_str(),s.size());		
		switch (member) {
			case FATHER:
				RG_father.insert(rgbuffer);
				break;
			case MOTHER:
				RG_mother.insert(rgbuffer);
				break;
			case SELF:
				RG_self.insert(rgbuffer);
				break;
			case SIBLING:
				RG_sibling.insert(rgbuffer);
				break;
			default:
				readgroups.insert(rgbuffer);
		}
	}
	
	if(readgroups.empty()) { // insert empty symbol
		readgroups.insert("null");
	}

	xfclose(fp);
}


// processGraph
//////////////////////////////////////////////////////////////////////////

void Microassembler::processGraph(Graph_t & g, const string & refname, const string & prefix, int minkmer, int maxkmer)
{	
	if (refname != "")
	{
		graphCnt++;
		
		//VERBOSE = false;
		
		// skip region if no mapped reads
		if(g.countMappedReads()<=0) { return; }

		cerr << "== Processing " << graphCnt << ": " << refname 
			<< " numsequences: " << g.readid2info.size() 
			<< " mapped: " << g.countMappedReads()
			<< " bastards: " << g.countBastardReads()
			<< endl;

		cerr << "=====================================================" << endl;

		// Load the reference
		map<string, Ref_t *>::iterator ri = reftable.find(refname);
		if (ri == reftable.end())
		{
			cerr << "Can't find ref info for: " << refname << endl;
			exit(1);
		}
		Ref_t * refinfo = ri->second;
		
		bool rptInRef = false;
		bool rptInQry = false;
		bool cycleInGraph = false;

		// dinamic kmer mode
		for (int k=minkmer; k<=maxkmer; k++) {
			g.setK(k);
			refinfo->setK(k);
			
			rptInRef = false;
			rptInQry = false;
			cycleInGraph = false;

			// exit if the region has a repeat of size K
			if(isRepeat(refinfo->rawseq, k)) { 
				cerr << "Repeat in reference sequence for kmer " << k << endl;
				rptInRef = true;
				//return; 
				continue;
			} 

			// exit if the region has an almost perfect repeat of size K
			if(isAlmostRepeat(refinfo->rawseq, k, MAX_MISMATCH)) { 
				cerr << "Near-perfect repeat in reference sequence for kmer " << k << endl;
				rptInRef = true;
				//return; 
				continue;
			}
			
			//if no repeats in the reference build graph			
			g.buildgraph(refinfo);
			
			double avgcov = ((double) g.totalreadbp_m) / ((double)refinfo->rawseq.length());			
			cerr << "reads: "   << g.readid2info.size()
				<< " reflen: "  << refinfo->rawseq.length()
				<< " readlen: " << g.totalreadbp_m
				<< " cov: "     << avgcov << endl;

			//printReads();
			g.printStats(0);
	
			string out_prefix = prefix + "/" + refname;
	
			if (PRINT_ALL) { g.printDot(out_prefix + ".0.dot",0); }
			
			// remove low covergae nodes and compute number of connected components
			g.removeLowCov(false, 0);
			int numcomp = g.markRefNodes();
			//cerr << "Num components = " << numcomp << endl;
			
			// process each connected components
			for (int c=1; c<=numcomp; c++) { 
				
				char comp[21]; // enough to hold all numbers up to 64-bits
				sprintf(comp, "%d", c);
				
				g.printStats(c); 
				
				// mark source and sink
				g.markRefEnds(refinfo, c);
			
				if (PRINT_ALL) { g.printDot(out_prefix + ".1l.c" + comp + ".dot", c); }
			
				// skip this component (and go to next one) if no tumor specific kmer found
				if ( !(g.hasTumorOnlyKmer()) ) { continue; }
					
				// if there is a cycle in the graph skip analysis
				if (g.hasCycle()) { g.clear(false); cycleInGraph = true; break; }

				g.checkReadStarts(c);
	
				// Initial compression
				g.compress(c); 
				g.printStats(c);
				if (PRINT_ALL) { g.printDot(out_prefix + ".2c.c" + comp + ".dot",c); }

				// Remove low coverage
				g.removeLowCov(true, c);
				if (PRINT_ALL) { g.printDot(out_prefix + ".3l.c" + comp + ".dot",c); }

				// Remove tips
				g.removeTips(c);
				if (PRINT_ALL) { g.printDot(out_prefix + ".4t.c" + comp + ".dot",c); }

				// skip analysis if there is a cycle in the graph 
				if (g.hasCycle()) { g.clear(false); cycleInGraph = true; break; }

				// skip analysis if there is a perfect or near-perfect repeat in the graph paths			
				if(g.hasRepeatsInGraphPaths()) { g.clear(false); rptInQry = true; break; }
			
				// Thread reads
				// BUG: threding is off because creates problems if the the bubble is not covered (end-to-end) 
				// by the reads. This is particularly problematic for detecting denovo events
				//g.threadReads();
				//if (PRINT_ALL) { g.printDot(out_prefix + ".4thread.dot"); }
		
				// scaffold contigs
				if (SCAFFOLD_CONTIGS)
				{
					g.scaffoldContigs();
				}

				if (PRINT_DENOVO)
				{
					g.denovoNodes(out_prefix + ".denovo.fa", refname);
				}  

				if (PRINT_REFPATH)
				{
					//g.markRefNodes();
					g.countRefPath(out_prefix + ".paths.fa", refname, false);
					//g.printFasta(prefix + "." + refname + ".nodes.fa");
				}

				if (PRINT_ALL) { g.printDot(out_prefix + ".final.c" + comp + ".dot",c); }					
			}
			
			if (rptInQry || cycleInGraph) { continue; }
			
			break; // break loop if graph has been processed correctly
		}
		
		// clear graph at the end.
		g.clear(true);
		
		if(rptInRef) { cerr << " Found repeat in reference" << endl; }
		if(rptInQry) { cerr << " Found repeat in assembly" << endl; }	
		if(cycleInGraph) { cerr << " Found cycle in assembly" << endl; }	
		
		cerr << "FINISHED" << endl;
	}
}

// fastqAsm
//////////////////////////////////////////////////////////////////////////

void Microassembler::fastqAsm(Graph_t & g, const string & prefix)
{
	cerr << endl << endl;
	cerr << "fastqAsm" << endl;
	cerr << "=====================================================" << endl;

	string refname = "fastq";
	string out_prefix = prefix + "." + refname;

	g.printStats();
		
	if (PRINT_ALL) 
    {
      g.printDot(out_prefix + ".0c.dot", 0);
      g.printFasta(out_prefix + ".0c.fa");
    }

	if (PRINT_ALL) 
    { 
      g.printDot(out_prefix + ".0.dot", 0); 
      g.printFasta(out_prefix + ".0.fa"); 
    }
    
	// Initial compression
	g.compress(0); 
	g.printStats();

	if (PRINT_RAW) 
    {
      g.printDot(out_prefix + ".1c.dot", 0);
      g.printFasta(out_prefix + ".1c.fa");
    }

	// Remove low coverage
	g.removeLowCov(true,0);
	if (PRINT_ALL) 
    { 
      g.printDot(out_prefix + ".2l.dot", 0); 
      g.printFasta(out_prefix + ".2l.fa"); 
    }

	// Remove tips
	g.removeTips(0);
	if (PRINT_ALL) 
    { 
      g.printDot(out_prefix + ".3t.dot", 0); 
      g.printFasta(out_prefix + ".3t.fa"); 
    }

	// Thread reads
	/*
	g.threadReads();
	if (PRINT_ALL) 
    { 
      g.printDot(out_prefix + ".4thread.dot"); 
      g.printFasta(out_prefix + ".4thread.fa"); 
    }
	*/
	// scaffold contigs
	if (SCAFFOLD_CONTIGS)
	{
		g.scaffoldContigs();
	}

	g.printDot(out_prefix + ".final.dot", 0);
	g.printFasta(out_prefix + ".final.fa");
	g.printPairs(out_prefix + ".pairs.fa");
	g.clear(true);

	cerr << endl;
}

int Microassembler::run(int argc, char** argv)
{
	stringstream HEADER;
	HEADER << 
		"\nProgram: Lancet (micro-assembly somatic variant caller)\n"
		"Version: 0.1.1 (beta), October 16 2015\n"
		"Contact: Giuseppe Narzisi <gnarzisi@nygenome.org>\n";
	
	string USAGE = "\nUsage: Lancet [options] --tumor <BAM file> --normal <BAM file> --ref <FASTA file>\n [-h for full list of commands]\n\n";

	if (argc == 1)
	{
		cerr << HEADER.str() << USAGE;
		exit(0);
	}

	cerr.setf(ios::fixed,ios::floatfield);
	cerr.precision(1);

	stringstream helptext;
	helptext << HEADER.str() << USAGE <<
		"Required\n"
		"   --tumor, -t          <BAM file>    : BAM file of mapped reads for tumor\n"
		"   --normal, -n         <BAM file>    : BAM file of mapped reads for normal\n"
		"   --ref, -r            <FASTA file>  : multifasta file of reference regions\n"
		"   --prefix, -p         <string>      : use prefix (default: mapfile)\n"
		"\nOptional\n"
		"   --min-k, k           <int>         : min kmersize (default: " << minK << ")\n"
		"   --max-k, -K          <int>         : max kmersize (default: " << maxK << ")\n"
		"   --trim-lowqual, -q   <int>         : trim bases below qv at 5' and 3' (default: " << MIN_QV << ")\n"
		"   --quality-range, -Q  <char>        : quality value range (default: " << (char) QV_RANGE << ")\n"
		"   --min-map-qual, -b   <inr>         : minimum read mapping quality in Phred-scale (default: " << MIN_MAP_QUAL << ")\n"
		"   --tip-len, -l        <int>         : max tip length (default: " << MAX_TIP_LEN << ")\n"
		"   --cov-thr, -c        <int>         : coverage threshold (default: " << COV_THRESHOLD << ")\n"
		"   --cov-ratio, -x      <float>       : minimum coverage ratio (default: " << MIN_COV_RATIO << ")\n"
		"   --max-avg-cov, -u    <int>         : maximum average coverage allowed per region (default: " << MAX_AVG_COV << ")\n"
		"   --tip-cov, -d        <int>         : tip coverage threshold (default: " << TIP_COV_THRESHOLD << ")\n"
		"   --dfs-limit, -F      <int>         : limit dfs search space (default: " << DFS_LIMIT << ")\n"
		"   --path-limit, -P     <int>         : limit on number of paths to report (default: " << PATH_LIMIT << ")\n"
		"   --max-indel-len, -T  <int>         : limit on size of detectable indel (default: " << MAX_INDEL_LEN << ")\n"
		"   --max-mismatch, -M   <int>         : max number of mismatches for near-perfect repeats (default: " << MAX_MISMATCH << ")\n"
		"   --rg-file, -g        <string>      : read group file\n"
		
		"\nFlags\n"
		"   -C            : extract reads from a BamFile \n"
		"   -D            : print de novo mutations (read map must be sorted and tagged with id)\n"
		"   -R            : print reference paths\n"
		"   -A            : print graph after every stage\n"
		"   -I            : don't print initial graph\n"
		"   -B            : include bastards\n"	
		"   -L <len>      : length of sequence to display at graph node (default: " << NODE_STRLEN << ")\n"
		"\n"
		"   -E            : fastq assembly, map file is in fq (experimental)\n"
		"   -v            : be verbose\n"
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
		{"prefix",  required_argument, 0, 'p'},
		{"rg-file",  required_argument, 0, 'g'},
		{"min-k",    required_argument, 0, 'k'},
		{"max-k",    required_argument, 0, 'K'},
		{"tip-len",  required_argument, 0, 'l'},
		{"cov-thr",  required_argument, 0, 'c'},
		{"cov-ratio",  required_argument, 0, 'x'},
		{"tip-cov",  required_argument, 0, 'd'},
		{"max-avg-cov",  required_argument, 0, 'u'},
		{"min-map-qual",  required_argument, 0, 'b'},
		{"trim-lowqual",  required_argument, 0, 'q'},
		{"quality-range",  required_argument, 0, 'Q'},

		{"node-str-len",  required_argument, 0, 'L'},
		{"dfs-limit",  required_argument, 0, 'F'},
		{"path-limit",  required_argument, 0, 'P'},
		{"max-indel-len",  required_argument, 0, 'T'},
		{"max-mismatch",  required_argument, 0, 'M'},

		{"erroflag", no_argument,      0, 'h'},		
		{"verbose", no_argument,       0, 'v'},
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
	while (!errflg && ((ch = getopt_long (argc, argv, "u:m:n:r:g:k:K:l:t:c:d:x:BDRACIhSL:T:M:vF:q:b:Q:P:p:E", long_options, &option_index)) != -1))
	{
		switch (ch)
		{
			case 't': TUMOR            = optarg;       break; 
			case 'n': NORMAL           = optarg;       break; 
			case 'r': REFFILE          = optarg;       break;
			
			case 'p': PREFIX           = optarg;       break;
			case 'g': RG_FILE          = optarg;       break;
			case 'k': minK             = atoi(optarg); break;
			case 'K': maxK             = atoi(optarg); break;
			case 'l': MAX_TIP_LEN      = atoi(optarg); break;
			//case 't': MIN_THREAD_READS = atoi(optarg); break;
			case 'c': COV_THRESHOLD    = atoi(optarg); break;
			case 'x': MIN_COV_RATIO    = atof(optarg); break;
			case 'd': TIP_COV_THRESHOLD= atoi(optarg); break;
			case 'u': MAX_AVG_COV      = atoi(optarg); break;
			
			case 'q': MIN_QV           = atoi(optarg); break;
			case 'b': MIN_MAP_QUAL     = atoi(optarg); break;
			case 'Q': QV_RANGE         = *optarg;      break;

			case 'L': NODE_STRLEN      = atoi(optarg); break;
			case 'F': DFS_LIMIT        = atoi(optarg); break;
			case 'P': PATH_LIMIT       = atoi(optarg); break;
			case 'T': MAX_INDEL_LEN    = atoi(optarg); break;
			case 'M': MAX_MISMATCH     = atoi(optarg); break;

			case 'C': BAMFILE          = 1; 		   break;
			case 'E': FASTQ_ASM        = 1;            break;
			case 'B': INCLUDE_BASTARDS = 1;            break;
			case 'v': VERBOSE          = 1;            break;
			case 'D': PRINT_DENOVO     = 1;            break;
			case 'R': PRINT_REFPATH    = 1;            break;
			case 'A': PRINT_ALL        = 1;            break;
			case 'I': PRINT_RAW        = 0;            break;

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

	if (FASTQ_ASM)
    {
		if (TUMOR == "") { cerr << "ERROR: Must provide a fastq file (-m)" << endl; errflg++; }
        if (PREFIX.empty()) { PREFIX = TUMOR; }
    }
	else
	{
		if (TUMOR == "") { cerr << "ERROR: Must provide the tumor BAM file (-m)" << endl; errflg++; }
		if (NORMAL == "") { cerr << "ERROR: Must provide the normal BAM file (-m)" << endl; errflg++; }
		if (REFFILE == "") { cerr << "ERROR: Must provide a reffile (-r)" << endl; errflg++; }
	}

	if (errflg) { exit(EXIT_FAILURE); }

	printConfiguration(cerr);

	if (REFFILE != "")
	{
		loadRefs(REFFILE);
	}

	// Process the reads
	FILE * fp = NULL;
	
	BamReader readerT;
	SamHeader headerT;
	RefVector referencesT;

	BamReader readerN;	
	SamHeader headerN;
	RefVector referencesN;
	
	
	if(BAMFILE) { // attempt to open our BamMultiReader
			if ( !readerT.Open(TUMOR) ) {
			cerr << "Could not open tumor BAM files." << endl;
			return -1;
		}
		// retrieve 'metadata' from BAM files, these are required by BamWriter
		headerT = readerT.GetHeader();
		referencesT = readerT.GetReferenceData();
		readerT.LocateIndex(); // locate and load BAM index file

		if ( !readerN.Open(NORMAL) ) {
			cerr << "Could not open normal BAM files." << endl;
			return -1;
		}
		// retrieve 'metadata' from BAM files, these are required by BamWriter
		headerN = readerN.GetHeader();
		referencesN = readerN.GetReferenceData();
		readerN.LocateIndex(); // locate and load BAM index file
	
		/*
		string outputBamFilename = PREFIX + "/outputBam.bam"; 
		// attempt to open our BamWriter
		if ( !writer.Open(outputBamFilename, header, references) ) {
			cerr << "Could not open output BAM file" << endl;
			return -1;
		}
		*/
		
		/*
		vector<RefData>::iterator it;
		cout << "refvector contains:" << endl;
		for ( it=references.begin() ; it < references.end(); it++ ) {
		cout << (*it).RefName << " " << (*it).RefLength << endl;
		}
		*/

		//load the read group information
		if(RG_FILE.compare("") != 0) {
			loadRG(RG_FILE,0);
		}
		else {
			readgroups.insert("null");
		}
		
		//set<string>::iterator prova;
		//for ( prova=readgroups.begin() ; prova != readgroups.end(); prova++ ) {
		//	cout << (*prova) <<  endl;
		//}
	}
	else { // text file with mapped reads
		fp = xfopen(TUMOR.c_str(), "r");
	}

	Graph_t g;

	//set configuration parameters
	g.setK(minK);
	g.setVerbose(VERBOSE);
	g.setMinQual(MIN_QUAL);
	g.setIncludeBastards(INCLUDE_BASTARDS);
	g.setBufferSize(BUFFER_SIZE);
	g.setDFSLimit(DFS_LIMIT);
	g.setPathLimit(PATH_LIMIT);
	g.setCovThreshold(COV_THRESHOLD);
	g.setMinCovRatio(MIN_COV_RATIO);
	g.setTipCovThreshold(TIP_COV_THRESHOLD);
	g.setPrintDotReads(PRINT_DOT_READS);
	g.setNodeStrlen(NODE_STRLEN);
	g.setMaxTipLength(MAX_TIP_LEN);
	g.setMaxIndelLen(MAX_INDEL_LEN);
	g.setMinThreadReads(MIN_THREAD_READS);
	g.setScaffoldContigs(SCAFFOLD_CONTIGS);
	g.setInsertSize(INSERT_SIZE);
	g.setInsertStdev(INSERT_STDEV);
	g.setMaxMismatch(MAX_MISMATCH);

	string graphref = "";

	char set  [BUFFER_SIZE];
	char code [BUFFER_SIZE];
	char chr  [BUFFER_SIZE];
	int  refstart;
	int  refend;
	char readname [BUFFER_SIZE];
	char seq1 [BUFFER_SIZE];
	char qv1  [BUFFER_SIZE];
	char seq2 [BUFFER_SIZE];
	char qv2  [BUFFER_SIZE];

	char refbuf[BUFFER_SIZE];

	int paircnt = 0;
	int graphcnt = 0;
	int readcnt = 0;

	int idx;

	if (FASTQ_ASM)
	{
		while (fscanf(fp, "%s", readname) == 1)
		{
			fscanf(fp, "%s", seq1);
			fscanf(fp, "%s", seq2); // header2
			fscanf(fp, "%s", qv1);

			g.addUnpaired(READSET, readname+1, seq1, qv1, Graph_t::CODE_MAPPED, UNDEFINED);
			readcnt++;
		}

		if (REFFILE == "")
		  {
		    fastqAsm(g, PREFIX);
		  }
		else
		  {
		    graphref = reftable.begin()->first;
		    processGraph(g, graphref, PREFIX, minK, maxK);
		  }

		graphcnt++;
	}
	else if (PRINT_DENOVO)
	{
		// To find denovos, the reads will be tagged with read set and sorted together
		while (fscanf(fp, "%s\t%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s",
			set, code, &idx, chr, &refstart, &refend, readname, seq1, qv1, seq2, qv2) == 11)
		{
			paircnt++;
			readcnt+=2;

			// printf("c:%s s:%d e:%d n:%s s:%s q:%s s:%s q:%s\n", 
			//         chr, refstart, refend, readname, seq1, qv1, seq2, qv2);
			//
			snprintf(refbuf, BUFFER_SIZE, "%s:%d-%d", chr, refstart, refend);

			string refstr = refbuf;

			if (refstr != graphref)
			{
				processGraph(g, graphref, PREFIX, minK, maxK);
				graphref = refstr;
				graphcnt++;
			}

			g.addPair(set, readname, seq1, qv1, seq2, qv2, code[0], UNDEFINED);
		}

		processGraph(g, graphref, PREFIX, minK, maxK);
	}
	else if(BAMFILE) {
		
		// for each reference location
		BamRegion region;
		map<string, Ref_t *>::iterator ri;
		for ( ri=reftable.begin() ; ri != reftable.end(); ri++ ) {
			graphref = (*ri).first;
			//cout << graphref << endl;
			Ref_t * refinfo = (*ri).second;
			
			// continue if the region has only Ns or prefect repeat of size maxK
			if(isNseq(refinfo->rawseq)) { continue; } 
			if(isRepeat(refinfo->rawseq, maxK)) { continue; } 

			region.LeftRefID = readerT.GetReferenceID(refinfo->refchr); // atoi((refinfo->refchr).c_str());
			region.RightRefID = readerT.GetReferenceID(refinfo->refchr); // atoi((refinfo->refchr).c_str());
			region.LeftPosition = refinfo->refstart;
			region.RightPosition = refinfo->refend;
			//cout << "region = " << refinfo->refchr << ":" << refinfo->refstart << "-" << refinfo->refend << endl; 

			bool jumpT = readerT.SetRegion(region);
			if(!jumpT) {
				cout << "Error: not able to jump successfully to the region's left boundary in tumor" << endl;
				return -1;
			}

			bool jumpN = readerN.SetRegion(region);
			if(!jumpN) {
				cout << "Error: not able to jump successfully to the region's left boundary in normal" << endl;
				return -1;
			}
			
			// iterate through all alignments
			//int num_PCR_duplicates = 0;
			BamAlignment al;
			string rg = "";
			int num_unmapped = 0;
			int totalreadbp = 0;
			double avgcov = 0.0;
			bool skip = false;
			
			
			/*** TUMOR ****/
			//while ( reader.GetNextAlignment(al) ) { // get next alignment and populate the alignment's string data fields
			while ( readerT.GetNextAlignmentCore(al) ) { // get next alignment and populate the alignment's string data fields
				
				avgcov = ((double) totalreadbp) / ((double)refinfo->rawseq.length());
				if(avgcov > MAX_AVG_COV) { 
					cout << "WARINING: Skip region " << refinfo->refchr << ":" << refinfo->refstart << "-" << refinfo->refend << ". Too much coverage (>" << MAX_AVG_COV << "x)." << endl;
					skip = true;
					break;
				}
				
				// skip alignments outside region
				int alstart = al.Position;
				int alend = al.GetEndPosition();
				if( (alstart < region.LeftPosition) || (alend > region.RightPosition) ) { continue; }
				
				if ( (al.MapQuality >= MIN_MAP_QUAL) && !al.IsDuplicate() ) { // only keeping ones with high map quality and skip PCR duplicates
					
					al.BuildCharData(); // Populates alignment string fields (read name, bases, qualities, tag data)
										
					int mate = 0;
					if(al.IsFirstMate()) { mate = 1; }
					if(al.IsSecondMate()) { mate = 2; }
				
					al.GetTag("RG", rg); // get the read group information for the read
					if(rg.empty()) { rg = "null"; }
					
					if ( (readgroups.find("null") != readgroups.end())  || (readgroups.find(rg) != readgroups.end()) ) { // select reads in the read group RG
						
						//writer.SaveAlignment(al); // save alignment to output bam file
						
						if (mate>1) { // mated pair
							if( !(al.IsMapped()) ) { // unmapped read
								g.addpaired("tumor", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_BASTARD, TMR);
								num_unmapped++; 
							}
							else { // mapped reads
								g.addpaired("tumor", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_MAPPED, TMR);								
							}
						}
						else { // unpaired
							g.addUnpaired("tumor", al.Name, al.QueryBases, al.Qualities, Graph_t::CODE_MAPPED, TMR);	
						}
						//cout << al.Name << endl;
						readcnt++;
						totalreadbp += (al.QueryBases).length();
						
						//void addMates(ReadId_t r1, ReadId_t r2)
						//{
						//g.readid2info[r1].mateid_m = r2;
						//g.readid2info[r2].mateid_m = r1;
						//}
					}
				}
				//else{ num_PCR_duplicates++; }
			}
			
			
			/*** NORMAL ****/
			//while ( reader.GetNextAlignment(al) ) { // get next alignment and populate the alignment's string data fields
			while ( readerN.GetNextAlignmentCore(al) ) { // get next alignment and populate the alignment's string data fields
				
				avgcov = ((double) totalreadbp) / ((double)refinfo->rawseq.length());
				if(avgcov > MAX_AVG_COV) { 
					cout << "WARINING: Skip region " << refinfo->refchr << ":" << refinfo->refstart << "-" << refinfo->refend << ". Too much coverage (>" << MAX_AVG_COV << "x)." << endl;
					skip = true;
					break;
				}
				
				// skip alignments outside region
				int alstart = al.Position;
				int alend = al.GetEndPosition();
				if( (alstart < region.LeftPosition) || (alend > region.RightPosition) ) { continue; }
				
				if ( (al.MapQuality >= MIN_MAP_QUAL) && !al.IsDuplicate() ) { // only keeping ones with high map quality and skip PCR duplicates
					
					al.BuildCharData(); // Populates alignment string fields (read name, bases, qualities, tag data)
										
					int mate = 0;
					if(al.IsFirstMate()) { mate = 1; }
					if(al.IsSecondMate()) { mate = 2; }
				
					al.GetTag("RG", rg); // get the read group information for the read
					if(rg.empty()) { rg = "null"; }
					
					if ( (readgroups.find("null") != readgroups.end())  || (readgroups.find(rg) != readgroups.end()) ) { // select reads in the read group RG
						
						//writer.SaveAlignment(al); // save alignment to output bam file
						
						if (mate>1) { // mated pair
							if( !(al.IsMapped()) ) { // unmapped read
								g.addpaired("normal", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_BASTARD, NML);
								num_unmapped++; 
							}
							else { // mapped reads
								g.addpaired("normal", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_MAPPED, NML);
							}
						}
						else { // unpaired
							g.addUnpaired("normal", al.Name, al.QueryBases, al.Qualities, Graph_t::CODE_MAPPED, NML);	
						}
						//cout << al.Name << endl;
						readcnt++;
						totalreadbp += (al.QueryBases).length();
						
						//void addMates(ReadId_t r1, ReadId_t r2)
						//{
						//g.readid2info[r1].mateid_m = r2;
						//g.readid2info[r2].mateid_m = r1;
						//}

					}
				}
				//else{ num_PCR_duplicates++; }
			}
			
			// close the reader & writer
			//cout << "Number of PCR duplicates: " << num_PCR_duplicates << endl;	
			//cout << "Number of unmapped reads: " << num_unmapped << endl;	
			
			if(!skip){ 
				processGraph(g, graphref, PREFIX, minK, maxK); 
			}
			else { g.clear(true); }
		}
		readerT.Close();
		readerN.Close();
	}
	else
	{
		// single set of reads against a reference

		while (fscanf(fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\t%s\t%s",
				code, &idx, chr, &refstart, &refend, readname, seq1, qv1, seq2, qv2) == 10)
		{
			paircnt++;
			readcnt+=2;

			// printf("c:%s s:%d e:%d n:%s s:%s q:%s s:%s q:%s\n", 
			//         chr, refstart, refend, readname, seq1, qv1, seq2, qv2);
			//
			snprintf(refbuf, BUFFER_SIZE, "%s:%d-%d", chr, refstart, refend);

			string refstr = refbuf;

			if (refstr != graphref)
			{
				processGraph(g, graphref, PREFIX, minK, maxK);
				graphref = refstr;
				graphcnt++;
			}
			g.addPair(READSET, readname, seq1, qv1, seq2, qv2, code[0], UNDEFINED);
		}
		processGraph(g, graphref, PREFIX, minK, maxK);
	}

	cerr << "=======" << endl;

	cerr << "total reads: " << readcnt << " pairs: " << paircnt << " total graphs: " << graphcnt << " ref sequences: " << reftable.size() <<  endl;

	if(!BAMFILE) {
		if (!feof(fp))
		{
			cerr << "ERROR: map file not at EOF." << endl;

			fgets(refbuf, BUFFER_SIZE, fp);
			fprintf(stderr, "line: %s\n", refbuf);

			fgets(refbuf, BUFFER_SIZE, fp);
			fprintf(stderr, "line: %s\n", refbuf);

			if (PRINT_DENOVO)
			{
				cerr << "WARNING: map file must have read sets tagged and sorted together for finding de novos" << endl;
			}
		}
	}
	return 0;
}

// main
//////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv)
{
	try {
		Microassembler* assembler = new Microassembler();
		assembler->run(argc, argv);
	}
	catch (int e) {
		cout << "An exception occurred. Exception Nr. " << e << endl;
	} 
	//catch(std::out_of_range& e) {
   	//	cerr << e.what( ) << endl;
 	//}
	//catch (...) { cout << "default exception"; }
}
