#include "Microassembler.hh"

/****************************************************************************
** Microassembler.cc
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

// load Red Groups
//////////////////////////////////////////////////////////////

void Microassembler::loadRG(const string & filename, int member) {

	FILE * fp;
	
	fp = xfopen(filename, "r");
	
	char rgbuffer[BUFFER_SIZE];
	string rg;

	while (fscanf(fp, "%s\n", rgbuffer) == 1) {
		//memcpy(a,s.c_str(),s.size());		
		readgroups.insert(rgbuffer);
	}
	
	if(readgroups.empty()) { // insert empty symbol
		readgroups.insert("null");
	}

	xfclose(fp);
}

// extract sammple name from SamHeader
//////////////////////////////////////////////////////////////
string Microassembler::retriveSampleName(SamHeader &header) {
	
	string sample_name = "NA";
	
	if(header.HasReadGroups()) {
		SamReadGroupDictionary RGS = header.ReadGroups;
		if(!(RGS.IsEmpty())) {
			SamReadGroup RG = *(RGS.Begin());
			if(RG.HasSample()) {
				sample_name = RG.Sample;
				//cerr << sample_name << endl;
			}
		}
	}
	return sample_name;
}


// processGraph
//////////////////////////////////////////////////////////////////////////

void Microassembler::processGraph(Graph_t & g, const string & refname, int minkmer, int maxkmer)
{	
	if (refname != "")
	{
		graphCnt++;
		
		//VERBOSE = false;
		
		// skip region if no mapped reads
		if(g.countMappedReads()<=0) { return; }

		if(verbose) { 
			cerr << "== Processing " << graphCnt << ": " << refname 
			<< " numsequences: " << g.readid2info.size() 
			<< " mapped: " << g.countMappedReads()
			<< " bastards: " << g.countBastardReads()
			<< endl;
			cerr << "=====================================================" << endl;
		}
		// Load the reference
		map<string, Ref_t *>::iterator ri = reftable->find(refname);
		if (ri == reftable->end())
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
				if(verbose) { cerr << "Repeat in reference sequence for kmer " << k << endl; }
				rptInRef = true;
				//return; 
				continue;
			} 

			// exit if the region has an almost perfect repeat of size K
			if(isAlmostRepeat(refinfo->rawseq, k, MAX_MISMATCH)) { 
				if(verbose) { cerr << "Near-perfect repeat in reference sequence for kmer " << k << endl; }
				rptInRef = true;
				//return; 
				continue;
			}
			
			//if no repeats in the reference build graph			
			g.buildgraph(refinfo);
			
			// error correct reads (just singletons)
			if(KMER_RECOVERY) {
				ErrorCorrector EC;
				EC.mersRecovery(g.nodes_m, 2, MIN_QUAL_CALL);
			}
			
			double avgcov = ((double) g.totalreadbp_m) / ((double)refinfo->rawseq.length());	
			if(verbose) {
				cerr << "reads: "   << g.readid2info.size()
				<< " reflen: "  << refinfo->rawseq.length()
				<< " readlen: " << g.totalreadbp_m
				<< " cov: "     << avgcov << endl;
			}
			//printReads();
			if(verbose) { g.printStats(0); }
	
			string out_prefix = "./" + refname;
	
			g.markRefNodes();
			if (PRINT_ALL) { g.printDot(out_prefix + ".0.dot",0); }
			
			// remove low covergae nodes and compute number of connected components
			g.removeLowCov(false, 0);
			int numcomp = g.markConnectedComponents();
			//cerr << "Num components = " << numcomp << endl;
			
			// process each connected components
			for (int c=1; c<=numcomp; c++) { 
				
				char comp[21]; // enough to hold all numbers up to 64-bits
				sprintf(comp, "%d", c);
				
				if(verbose) { g.printStats(c); }
				
				// mark source and sink
				g.markRefEnds(refinfo, c);
			
				if (PRINT_ALL) { g.printDot(out_prefix + ".1l.c" + comp + ".dot", c); }
			
				// skip this component (and go to next one) if no tumor specific kmer found
				//if ( !(g.hasTumorOnlyKmer()) ) { continue; }
					
				// if there is a cycle in the graph skip analysis
				if (g.hasCycle()) { g.clear(false); cycleInGraph = true; break; }

				g.checkReadStarts(c);
	
				// Initial compression
				g.compress(c); 
				if(verbose) { g.printStats(c); }
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
				if(g.hasRepeatsInGraphPaths(refinfo)) { g.clear(false); rptInQry = true; break; }
			
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

				/*
				if (PRINT_DENOVO)
				{
					g.denovoNodes(out_prefix + ".denovo.fa", refname);
				}  
				*/

				//g.markRefNodes();
				g.countRefPath(out_prefix + ".paths.fa", refname, false);
				//g.printFasta(prefix + "." + refname + ".nodes.fa");

				if (PRINT_ALL) { g.printDot(out_prefix + ".final.c" + comp + ".dot",c); }					
			}
			
			if (rptInQry || cycleInGraph) { continue; }
			
			break; // break loop if graph has been processed correctly
		}
		
		// clear graph at the end.
		g.clear(true);
		
		if(verbose) {
			if(rptInRef) { cerr << " Found repeat in reference" << endl; }
			if(rptInQry) { cerr << " Found repeat in assembly" << endl; }	
			if(cycleInGraph) { cerr << " Found cycle in assembly" << endl; }	
		}
		
		if(verbose) { cerr << "FINISHED" << endl; }
	}
}


// isActiveRegion
// Examines reads alignments (CIGAR and MD) to find evidence of mutations
// returns true if there is evidence of mutation in the region
bool Microassembler::isActiveRegion(BamReader &reader, Ref_t *refinfo, BamRegion &region, int code) {
	
	// iterate through all alignments
	int MIN_EVIDENCE = filters.minAltCntTumor; // min evidence equal to min support for somatic variant
	BamAlignment al;
	int totalreadbp = 0;
	bool ans = false;
	int MQ = MIN_MAP_QUAL;
	string CIGAR = "";
	string rg = "";
	string md = "";
	
	vector< int > clipSizes;
	vector< int > readPositions; 
	vector< int > genomePositions;

    map<int,int> mapX; // table with counts of all mismatches at a given locus	
    map<int,int> mapI; // table with counts of all insertions at a given locus	
    map<int,int> mapD; // table with counts of all deletion at a given locus	
    map<int,int> mapSC; // table with counts of all softclipped sequences starting at a given locus
    map<int,int>::iterator mit;
	
	// more sensitive in normal (extract all reads)
	if (code == NML) { MQ = 0; }
	
	/*** TUMOR ****/
	//while ( reader.GetNextAlignment(al) ) { // get next alignment and populate the alignment's string data fields
	while ( reader.GetNextAlignmentCore(al) ) { // get next alignment and populate the alignment's string data fields
		
		int alstart = al.Position;
		int alend = al.GetEndPosition();
		if( (alstart < region.LeftPosition) || (alend > region.RightPosition) ) { continue; } // skip alignments outside region
		
		if ( (al.MapQuality >= MQ) && !al.IsDuplicate() ) { // only keep reads with high map quality and skip PCR duplicates
			
			al.BuildCharData(); // Populates alignment string fields (read name, bases, qualities, tag data)
			
			al.GetTag("RG", rg); // get the read group information for the read
			if(rg.empty()) { rg = "null"; }
			
			if ( (readgroups.find("null") != readgroups.end())  || (readgroups.find(rg) != readgroups.end()) ) { // select reads in the read group RG
				
				// parse MD string
				// String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*10
				al.GetTag("MD", md); // get string of mismatching positions
				parseMD(md, mapX, alstart);
				
				// add SNV to database
				//Variant_t(string chr_, int pos_, string ref_, string alt_, int ref_cov_normal_, int ref_cov_tumor_, int alt_cov_normal_fwd_, int alt_cov_normal_rev_, int alt_cov_tumor_fwd_, int alt_cov_tumor_rev_, char prev_bp_ref_, char prev_bp_alt_, Filters &fs)
				
				
				// example: 31M1I17M1D37M
				CIGAR = "";
				int pos = alstart; // initialize position to start of alignment
				for (std::vector<CigarOp>::iterator it = (al.CigarData).begin() ; it != (al.CigarData).end(); ++it) {
					char T = (*it).Type;
					
					// update position (except for insertions)
					if(T!='I') { pos += (*it).Length; }
					
					if(T=='X') {
						//cerr << "X:" << pos << "|";
						mit = mapX.find(pos);
						if (mit != mapX.end()) { ((*mit).second)++; }
						else { mapX.insert(std::pair<int,int>(pos,1)); }
					}
					if(T=='I') {
						//cerr << "I:" << pos << "|";
						mit = mapI.find(pos);
						if (mit != mapI.end()) { ((*mit).second)++; }
						else { mapI.insert(std::pair<int,int>(pos,1)); }
					}
					if(T=='D') {
						//cerr << "D:" << pos << "|";
						mit = mapD.find(pos);
						if (mit != mapD.end()) { ((*mit).second)++; }
						else { mapD.insert(std::pair<int,int>(pos,1)); }
					}
					
					std::stringstream ss;
					ss << (*it).Length;
					CIGAR += ss.str();
					CIGAR += (*it).Type;
				}
				//cerr << endl;
				
				int len = al.Length;
				totalreadbp += len;
				
				clipSizes.clear();
				readPositions.clear();
				genomePositions.clear();
				
				//cerr << CIGAR << " MD:" << md << " START:" << alstart << " ";
				if(al.GetSoftClips(clipSizes, readPositions, genomePositions)) {	
					for (std::vector<int>::iterator it = genomePositions.begin() ; it != genomePositions.end(); ++it) {
						//cerr << (*it) << " ";
						mit = mapSC.find((*it));
						if (mit != mapSC.end()) { ((*mit).second)++; }
						else { mapSC.insert(std::pair<int,int>((*it),1)); }
					}
				}
				//cerr << endl;
			}
		}
	}

	// check for any locus with evidence for SNVs, indel, or soft-clipped sequences	
	bool snv_evidence = false;
	bool indel_evidence = false;
	bool softclip_evidence = false;
	
	//cerr << "X: ";
    for (mit=mapX.begin(); mit!=mapX.end(); ++mit) {
		if((*mit).second >= MIN_EVIDENCE) { 
			//cerr << (*mit).first << "," << (*mit).second << "|";
			snv_evidence = true;
			ans = true; 
			break; 
		}
	}
	//cerr << endl;
	
	//cerr << "I: ";
    for (mit=mapI.begin(); mit!=mapI.end(); ++mit) {
		if((*mit).second >= MIN_EVIDENCE) { 
			//cerr << (*mit).first << "," << (*mit).second << "|";
			indel_evidence = true;
			ans = true; 
			break; 
		}
	}
	//cerr << endl;
	
	//cerr << "D: ";
    for (mit=mapD.begin(); mit!=mapD.end(); ++mit) {
		if((*mit).second >= MIN_EVIDENCE) { 
			//cerr << (*mit).first << "," << (*mit).second << "|";
			indel_evidence = true;
			ans = true; 
			break;
		}
	}
	//cerr << endl;
	
	//cerr << "S: ";
    for (mit=mapSC.begin(); mit!=mapSC.end(); ++mit) {
		if((*mit).second >= MIN_EVIDENCE) { 
			//cerr << (*mit).first << "," << (*mit).second << "|";
			softclip_evidence = true;
			ans = true; 
			break; 
		}
	}
	//cerr << endl;
	
	if(code == TMR) {
		if(snv_evidence && !indel_evidence && !softclip_evidence)   { num_snv_only_regions++; ans = false; }
		if(!snv_evidence && indel_evidence && !softclip_evidence)   { num_indel_only_regions++; }
		if(!snv_evidence && !indel_evidence && softclip_evidence)   { num_softclip_only_regions++; }
		if(!snv_evidence && (indel_evidence || softclip_evidence))  { num_indel_or_softclip_regions++; }	
		if(snv_evidence || indel_evidence || softclip_evidence)     { num_snv_or_indel_or_softclip_regions++; }
	}

	if(snv_evidence && !indel_evidence && !softclip_evidence) { ans = false; }
	
	return ans;
}

// extract the reads from BAM file
// return false if the region could not be analyzed (e.g., too much coverage)
bool Microassembler::extractReads(BamReader &reader, Graph_t &g, Ref_t *refinfo, BamRegion &region, int &readcnt, int code) {
	
	// iterate through all alignments
	//int num_PCR_duplicates = 0;
	BamAlignment al;
	string rg = "";
	int num_unmapped = 0;
	int totalreadbp = 0;
	double avgcov = 0.0;
	int MQ = MIN_MAP_QUAL;
	bool skip = false;
	
	// more sensitive in normal (extract all reads)
	if (code == NML) { MQ = 0; }
		
	/*** TUMOR ****/
	//while ( reader.GetNextAlignment(al) ) { // get next alignment and populate the alignment's string data fields
	while ( reader.GetNextAlignmentCore(al) ) { // get next alignment and populate the alignment's string data fields
		
		avgcov = ((double) totalreadbp) / ((double)refinfo->rawseq.length());
		if(avgcov > MAX_AVG_COV) { 
			cerr << "WARNING: Skip region " << refinfo->refchr << ":" << refinfo->refstart << "-" << refinfo->refend << ". Too much coverage (>" << MAX_AVG_COV << "x)." << endl;
			skip = true;
			break;
		}
		
		int alstart = al.Position;
		int alend = al.GetEndPosition();
		if( (alstart < region.LeftPosition) || (alend > region.RightPosition) ) { continue; } // skip alignments outside region
		
		if ( (al.MapQuality >= MQ) && !al.IsDuplicate() ) { // only keep reads with high map quality and skip PCR duplicates
			
			al.BuildCharData(); // Populates alignment string fields (read name, bases, qualities, tag data)
								
			int mate = 0;
			int strand = FWD;
			if(al.IsFirstMate()) { mate = 1; }
			if(al.IsSecondMate()) { mate = 2; }
			if(al.IsReverseStrand()) { strand = REV; }
			
			/*
			string oq;
			al.GetTag("OQ", oq); // get original base quality scores
			cerr << oq << endl;
			*/
			
			al.GetTag("RG", rg); // get the read group information for the read
			if(rg.empty()) { rg = "null"; }
			
			if ( (readgroups.find("null") != readgroups.end())  || (readgroups.find(rg) != readgroups.end()) ) { // select reads in the read group RG
												
				if (mate>1) { // mated pair
					//cerr << "PAIRED!!" << endl;
					if( !(al.IsMapped()) ) { // unmapped read
						g.addpaired("tumor", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_BASTARD, code, strand);
						//g.addpaired("tumor", al.Name, al.QueryBases, oq, mate, Graph_t::CODE_BASTARD, code, strand);
						num_unmapped++; 
					}
					else { // mapped reads
						g.addpaired("tumor", al.Name, al.QueryBases, al.Qualities, mate, Graph_t::CODE_MAPPED, code, strand);								
						//g.addpaired("tumor", al.Name, al.QueryBases, oq, mate, Graph_t::CODE_MAPPED, code, strand);								
					}
				}
				else { // unpaired
					//cerr << "UNPAIRED!!" << endl;
					g.addUnpaired("tumor", al.Name, al.QueryBases, al.Qualities, Graph_t::CODE_MAPPED, code, strand);	
					//g.addUnpaired("tumor", al.Name, al.QueryBases, oq, Graph_t::CODE_MAPPED, code, strand);	
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
	
	return skip;
}

// extract the reads from BAMs and process them
int Microassembler::processReads() {
	
	cerr << "Process reads" << endl;
	
	// Process the reads
	BamReader readerT;
	SamHeader headerT;
	RefVector referencesT;

	BamReader readerN;	
	SamHeader headerN;
	RefVector referencesN;
			
	// attempt to open our BamMultiReader
	if ( !readerT.Open(TUMOR) ) {
		cerr << "Could not open tumor BAM files." << endl;
		return -1;
	}
	// retrieve 'metadata' from BAM files, these are required by BamWriter
	headerT = readerT.GetHeader();
	referencesT = readerT.GetReferenceData();
	readerT.LocateIndex(); // locate and load BAM index file
	sample_name_tumor = retriveSampleName(headerT); // extract sammple name

	if ( !readerN.Open(NORMAL) ) {
		cerr << "Could not open normal BAM files." << endl;
		return -1;
	}
	// retrieve 'metadata' from BAM files, these are required by BamWriter
	headerN = readerN.GetHeader();
	referencesN = readerN.GetReferenceData();
	readerN.LocateIndex(); // locate and load BAM index file
	sample_name_normal = retriveSampleName(headerN); // extract sammple name

	//load the read group information
	if(RG_FILE.compare("") != 0) {
		loadRG(RG_FILE,0);
	}
	else {
		readgroups.insert("null");
	}

	Graph_t g;

	//set configuration parameters
	g.setDB(&vDB);
	g.setK(minK);
	g.setVerbose(verbose);
	g.setMoreVerbose(VERBOSE);
	g.setMinQualTrim(MIN_QUAL_TRIM);
	g.setMinQualCall(MIN_QUAL_CALL);
	g.setBufferSize(BUFFER_SIZE);
	g.setDFSLimit(DFS_LIMIT);
	g.setCovThreshold(COV_THRESHOLD);
	g.setMinCovRatio(MIN_COV_RATIO);
	g.setLowCovThreshold(LOW_COV_THRESHOLD);
	g.setPrintDotReads(PRINT_DOT_READS);
	g.setNodeStrlen(NODE_STRLEN);
	g.setMaxTipLength(MAX_TIP_LEN);
	g.setMaxIndelLen(MAX_INDEL_LEN);
	g.setMinThreadReads(MIN_THREAD_READS);
	g.setScaffoldContigs(SCAFFOLD_CONTIGS);
	g.setInsertSize(INSERT_SIZE);
	g.setInsertStdev(INSERT_STDEV);
	g.setMaxMismatch(MAX_MISMATCH);
	g.setFilters(filters);

	string graphref = "";

	int paircnt = 0;
	int graphcnt = 0;
	int readcnt = 0;
	
	// for each reference location
	BamRegion region;
	map<string, Ref_t *>::iterator ri;
		
	int counter = 0;
	double progress;
	double old_progress = 0;
	for ( ri=reftable->begin() ; ri != reftable->end(); ri++ ) {

		counter++;
		progress = floor(100*(double(counter)/(double)reftable->size()));
		if (progress > old_progress) {
			cerr << "Thread " << ID << " is " << progress << "\% done." << endl;
			old_progress = progress;
		}
			
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
			cerr << "Error: not able to jump successfully to the region's left boundary in tumor" << endl;
			return -1;
		}

		bool jumpN = readerN.SetRegion(region);
		if(!jumpN) {
			cerr << "Error: not able to jump successfully to the region's left boundary in normal" << endl;
			return -1;
		}
		
		bool activeT = true;
		bool activeN = true;
		
		if (ACTIVE_REGION_MODULE) {
			activeT = isActiveRegion(readerT, refinfo, region, TMR);
			activeN = isActiveRegion(readerN, refinfo, region, NML);
		}
		
		if(activeT || activeN){
			
			readerT.SetRegion(region); // safe to jump back: errors would have been detected in the previous call to jump
			readerN.SetRegion(region); // safe to jump back: errors would have been detected in the previous call to jump
			
			bool skipT = extractReads(readerT, g, refinfo, region, readcnt, TMR);
			bool skipN = extractReads(readerN, g, refinfo, region, readcnt, NML);
			
			if(!skipT && !skipN) { 
				processGraph(g, graphref, minK, maxK);
			}
		}
		else {
			num_skip++;
			g.clear(true);
			if(verbose) { cerr << "Skip region: not enough evidence for variation." << endl; }
		}
	}
		
	readerT.Close();
	readerN.Close();
	
	if(verbose) cerr << "=======" << endl;
	if(verbose) cerr << "total reads: " << readcnt << " pairs: " << paircnt << " total graphs: " << graphcnt << " ref sequences: " << reftable->size() <<  endl;
	
	return 0;
}
