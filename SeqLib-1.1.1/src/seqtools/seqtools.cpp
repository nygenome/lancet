#include <getopt.h>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>

#include "SeqLib/BFC.h"
#include "SeqLib/FastqReader.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BWAWrapper.h"
#include "SeqLib/FermiAssembler.h"

void kt_pipeline(int n_threads, void *(*func)(void*, int, void*), void *shared_data, int n_steps);

#define AUTHOR "Jeremiah Wala <jwala@broadinstitute.org>"

static const char *SEQTOOLS_USAGE_MESSAGE =
"Program: seqtools \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools <command> [options]\n\n"
"Commands:\n"
"           bfc       Error correction from a BAM or fasta, direct to re-aligned BAM\n"
"           fml       FermiKit assembly (with error correction), direct to re-aligned BAM\n" 
"\nReport bugs to jwala@broadinstitute.org \n\n";

static const char *BFC_USAGE_MESSAGE =
"Program: seqtools bfc \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools bfc [options]\n\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --fasta,     -f        Output stream should be a FASTA (no realignment)\n"
"  --bam,       -b        Output stream should be a BAM (not SAM)\n"
"  --cram,      -C        Output stream should be a CRAM (not SAM)\n"
"  --infasta,   -F <file> Input a FASTA insted of BAM/SAM/CRAM stream\n"
"  --reference, -G <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

static const char *FML_USAGE_MESSAGE =
"Program: seqtools fml \n"
"Contact: Jeremiah Wala [ jwala@broadinstitute.org ]\n"
"Usage: seqtools fml [options]\n\n"
"Description: Extract sequences and assemble and realign contigs\n"
"Commands:\n"
"  --verbose,   -v        Set verbose output\n"
"  --fasta,     -f        Output stream should be a FASTA (no realignment)\n"
"  --bam,       -b        Output stream should be a BAM (not SAM)\n"
"  --cram,      -C        Output stream should be a CRAM (not SAM)\n"
"  --infasta,   -F <file> Input a FASTA insted of BAM/SAM/CRAM stream\n"
"  --reference, -G <file> Reference genome if using BWA-MEM realignment\n"
"\nReport bugs to jwala@broadinstitute.org \n\n";

void runbfc(int argc, char** argv);
void runfml(int argc, char** argv);
void parseOptions(int argc, char** argv, const char* msg);

namespace opt {

  static bool verbose = false;
  static char mode = 's';
  static std::string input;
  static std::string reference = "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta";
  static std::string fasta; // input is a fasta
  static std::string target; // input target sequence
}

static const char* shortopts = "hbfvCG:F:T:";
static const struct option longopts[] = {
  { "help",                    no_argument, NULL, 'h' },
  { "verbose",                 no_argument, NULL, 'v' },
  { "bam",                     no_argument, NULL, 'b' },
  { "cram",                    no_argument, NULL, 'C' },
  { "fasta",                   no_argument, NULL, 'f' },
  { "infasta",                 required_argument, NULL, 'F' },
  { "reference",               required_argument, NULL, 'G' },
  { "target",                  required_argument, NULL, 'T' },
  { NULL, 0, NULL, 0 }
};

int main(int argc, char** argv) {

   if (argc <= 1) {
    std::cerr << SEQTOOLS_USAGE_MESSAGE;
    return 0;
  } else {
    std::string command(argv[1]);
    if (command == "help" || command == "--help") {
      std::cerr << SEQTOOLS_USAGE_MESSAGE;
      return 0;
    } else if (command == "bfc") {
      runbfc(argc -1, argv + 1);
    } else if (command == "fml") {
      runfml(argc -1, argv + 1);
    } else {
      std::cerr << SEQTOOLS_USAGE_MESSAGE;
      return 0;
    }
  } 
  
  return 0;

}

void runfml(int argc, char** argv) {

  parseOptions(argc, argv, FML_USAGE_MESSAGE);
  
  SeqLib::FermiAssembler fml;

  if (!opt::fasta.empty()) {
    SeqLib::FastqReader f(opt::fasta);

    SeqLib::UnalignedSequenceVector usv;
    std::string qn, seq;
    SeqLib::UnalignedSequence u;
    if (opt::verbose)
      std::cerr << "...reading fasta/fastq file " << opt::fasta << std::endl;
    while (f.GetNextSequence(u)) 
      fml.AddRead(u);

    if (opt::verbose)
      std::cerr << "...read in " << SeqLib::AddCommas(fml.NumSequences()) << " sequences" <<  std::endl;
    
  } else {

    SeqLib::BamReader br;
    if (!br.Open(opt::input == "-" ? "-" : opt::input)) 
      exit(EXIT_FAILURE);
      
    if (opt::verbose)
      std::cerr << "...opened " << opt::input << std::endl;
    SeqLib::BamRecord rec;
    SeqLib::BamRecordVector brv;
    size_t count = 0;
    while(br.GetNextRecord(rec)) {
      if (++count % 1000000 == 0 && opt::verbose)
	std::cerr << "...at read " << SeqLib::AddCommas(count) << " " << rec.Brief() << std::endl;
      brv.push_back(rec); //rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
    }
    fml.AddReads(brv);

  }

  if (opt::verbose) 
    std::cerr << "...error correcting " << std::endl;
  fml.CorrectReads();

  if (opt::verbose) 
    std::cerr << "...assembling " << std::endl;
  fml.PerformAssembly();

  // retrieve the reads
  std::vector<std::string> contigs = fml.GetContigs();

  size_t count = 0;
  if (opt::mode == 'f') {
    for (std::vector<std::string>::const_iterator i = contigs.begin(); i != contigs.end(); ++i)
      std::cout << ">contig" << ++count << std::endl << *i << std::endl;
    return;
  }

  SeqLib::BamWriter bw;
  if (opt::mode == 'b')
    bw = SeqLib::BamWriter(SeqLib::BAM);
  else if (opt::mode == 's')
    bw = SeqLib::BamWriter(SeqLib::SAM);
  else if (opt::mode == 'C') {
    bw = SeqLib::BamWriter(SeqLib::CRAM);
    bw.SetCramReference(opt::reference);
  }
  else {
    std::cerr << "Unrecognized output stream mode " << opt::mode << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (!bw.Open("-"))
    exit(EXIT_FAILURE);
  
  SeqLib::BWAWrapper bwa;
  if (!bwa.LoadIndex(opt::reference)) {
    std::cerr << "Failed to load index for BWA-MEM from: " << opt::reference << std::endl;
    exit(EXIT_FAILURE);
  }
  
  bw.SetHeader(bwa.HeaderFromIndex());
  bw.WriteHeader();

  if (opt::verbose) 
    std::cerr << "...aligning contig with BWA-MEM" << std::endl;
 
  // run through and read
  std::stringstream ss;
  for (std::vector<std::string>::const_iterator i = contigs.begin(); i != contigs.end(); ++i) {
    SeqLib::BamRecordVector brv;
    const bool hardclip = false;
    const double frac = 0.9;
    ss << "contig" << ++count;
    const int max_secondary = 10;
    bwa.AlignSequence(*i, ss.str(), brv, hardclip, frac, max_secondary);
    ss.str(std::string());
    for (SeqLib::BamRecordVector::iterator r = brv.begin();
	 r != brv.end(); ++r) {
      bw.WriteRecord(*r);
    }
  }
  
}

void runbfc(int argc, char** argv) {

  parseOptions(argc, argv, BFC_USAGE_MESSAGE);

  SeqLib::BFC b;

  if (!opt::fasta.empty()) {
    // read in a fasta file
    SeqLib::FastqReader f(opt::fasta);
    
    SeqLib::UnalignedSequence u;
    while (f.GetNextSequence(u)) {
      if (!b.AddSequence(u.Seq.c_str(), u.Qual.c_str(), u.Name.c_str())) {
	std::cerr << "Error adding sequence from fasta: " << u.Seq << std::endl;
	exit(EXIT_FAILURE);
      }
    }
  } else { //if (opt::mode == 'b' || opt::mode == 's' || opt::mode == 'C') {
    SeqLib::BamReader br;
    if (!br.Open(opt::input == "-" ? "-" : opt::input))
      exit(EXIT_FAILURE);
    if (opt::verbose)
      std::cerr << "...opened " << opt::input << std::endl;
    SeqLib::BamRecord rec;
    size_t count = 0;
    while(br.GetNextRecord(rec)) {
      if (++count % 1000000 == 0 && opt::verbose)
	std::cerr << "...at read " << SeqLib::AddCommas(count) << " " << rec.Brief() << std::endl;
      //b.AddSequence(rec); //rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
      b.AddSequence(rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
    }
  } 

  if (opt::verbose)
    std::cerr << "...read in " << SeqLib::AddCommas(b.NumSequences()) << " sequences" << std::endl;
  
  if (!b.Train()) {
    std::cerr << "Training failed on " << b.NumSequences() << std::endl;
    exit(EXIT_FAILURE);
  }
  if (opt::verbose)
    std::cerr << "...finished training " << SeqLib::AddCommas(b.NumSequences()) << " sequences" << std::endl;
  if (!b.ErrorCorrect()) {
    std::cerr << "Correction failed on " << b.NumSequences() << std::endl;
    exit(EXIT_FAILURE);
  }
  if (opt::verbose)
    std::cerr << "...finished correcting " << SeqLib::AddCommas(b.NumSequences()) << " sequences" << std::endl;

  SeqLib::UnalignedSequenceVector u;
  std::string seqr, namr;
  while (b.GetSequence(seqr, namr)) 
    u.push_back(SeqLib::UnalignedSequence(namr, seqr));

  if (opt::verbose)
    std::cerr << "nseqs: " << u.size() 
	      << " kcov: " << b.GetKCov() 
	      << " kmer: " << b.GetKMer() << std::endl;
  
  if (opt::mode == 'f') {
    for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin();
	 i != u.end(); ++i) {
      std::cout << ">" << i->Name << std::endl << i->Seq << std::endl;
    }
    return;
  } 

  SeqLib::BamWriter bw;
  if (opt::mode == 'b')
    bw = SeqLib::BamWriter(SeqLib::BAM);
  else if (opt::mode == 's')
    bw = SeqLib::BamWriter(SeqLib::SAM);
  else if (opt::mode == 'C') {
    bw = SeqLib::BamWriter(SeqLib::CRAM);
    bw.SetCramReference(opt::reference);
  }
  else {
    std::cerr << "Unrecognized output stream mode " << opt::mode << std::endl;
    exit(EXIT_FAILURE);
  }
  
  if (!bw.Open("-"))
    exit(EXIT_FAILURE);
  
  SeqLib::BWAWrapper bwa;
  if (opt::verbose)
    std::cerr << "...loading reference genome" << std::endl;
  if (!bwa.LoadIndex(opt::reference)) {
    std::cerr << "Failed to load index for BWA-MEM from: " << opt::reference << std::endl;
    exit(EXIT_FAILURE);
  }
  
  bw.SetHeader(bwa.HeaderFromIndex());
  bw.WriteHeader();

  if (opt::verbose)
    std::cerr << "...realigning corrected sequences with BWA-MEM" << std::endl;
  // run through and read
  for (SeqLib::UnalignedSequenceVector::const_iterator i = u.begin(); i != u.end(); ++i) {
    SeqLib::BamRecordVector brv;
    const bool hardclip = false;
    const double frac = 0.9;
    const int max_secondary = 10;
    bwa.AlignSequence(i->Seq, i->Name, brv, hardclip, frac, max_secondary);
    for (SeqLib::BamRecordVector::iterator r = brv.begin();
	 r != brv.end(); ++r) {
      if (!i->Qual.empty())
	r->SetQualities(i->Qual, 33);
      bw.WriteRecord(*r);
    }
  }
  
}
// parse the command line options
void parseOptions(int argc, char** argv, const char* msg) {

  bool die = false;
  bool help = false;

  // get the first argument as input
  if (argc > 1)
    opt::input = std::string(argv[1]);
  
  for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'v': opt::verbose = true; break;
    case 'f': opt::mode = 'f'; break;
    case 'F': arg >> opt::fasta; break;
    case 'b': opt::mode = 'b'; break;
    case 'C': opt::mode = 'C'; break;
    case 'T': arg >> opt::target; break;
    case 'G': arg >> opt::reference; break;
    default: die= true; 
    }
  }

  if (die || help || (opt::input.empty() && opt::fasta.empty())) {
      std::cerr << "\n" << msg;
      if (die)
	exit(EXIT_FAILURE);
      else 
	exit(EXIT_SUCCESS);	
    }
}
