#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include<boost/test/unit_test.hpp>

#include <climits>
#include <boost/test/unit_test.hpp>

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/ReadFilter.h"
#include "SeqLib/FermiAssembler.h"
#include "SeqLib/SeqPlot.h"
#include "SeqLib/RefGenome.h"

#define GZBED "test_data/test.bed.gz"
#define GZVCF "test_data/test.vcf.gz"
#define SBAM "test_data/small.bam"
#define OBAM "test_data/small_out.bam"
#define OCRAM "test_data/small_out.cram"
#define HGREF "/seq/references/Homo_sapiens_assembly19/v1/Homo_sapiens_assembly19.fasta"
#define TREF "test_data/test_ref.fa"
#define OREF "tmp_output.fa"
#define BEDFILE "test_data/test.bed"
#define VCFFILE "test_data/test.vcf"
#define JSON1 "test_data/example4.json"

using namespace SeqLib::Filter;
using namespace SeqLib;

#include <fstream>
#include "SeqLib/BFC.h"

BOOST_AUTO_TEST_CASE( read_gzbed ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");

  SeqLib::GRC g(GZBED, br.Header());

  BOOST_CHECK_EQUAL(g.size(), 3);

  BOOST_CHECK_EQUAL(g[2].chr, 21);

  SeqLib::GRC v(GZVCF, br.Header());
  BOOST_CHECK_EQUAL(v.size(), 57);

  BOOST_CHECK_EQUAL(v[29].chr, 0);
}

BOOST_AUTO_TEST_CASE ( bfc ) {

  BFC b;

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");

  SeqLib::BamRecord rec;
  BamRecordVector brv, brv2;
  size_t count = 0;
  while(br.GetNextRecord(rec) && count++ < 10000) 
    brv.push_back(rec);

  for (int i = 5000; i < brv.size(); ++i)
    brv2.push_back(brv[i]);

  count = 0;
  while(br.GetNextRecord(rec) && count++ < 10000) 
    brv2.push_back(rec);

  std::ofstream orig("orig.fa");
  std::ofstream corr("corr.fa");
    
  for (auto& i : brv)
    orig << ">" << i.Qname() << std::endl << i.Sequence() << std::endl;

  // add the seqs
  for (auto& r : brv)
    b.AddSequence(r.Sequence().c_str(), r.Qualities().c_str(), r.Qname().c_str());

  b.Train();
  b.clear();

  //b.ErrorCorrectToTag(brv2, "KC");  

  //UnalignedSequenceVector v;
  //b.GetSequences(v);

  // write to corrected
  //for (auto& i : v) {
  //  corr << ">" << i.Name << std::endl << i.Seq << std::endl;
  //}
  //orig.close();
  //corr.close();
  
  //
  //v.clear();
  //b.FilterUnique();
  //b.GetSequences(v);

  // do everything at once
  //b.TrainAndCorrect(brv2);

  // do everything in place
  //b.TrainCorrection(brv2);
  //b.ErrorCorrectInPlace(brv2);
}

BOOST_AUTO_TEST_CASE( correct_and_assemble ) {

  BFC b;

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");

  SeqLib::BamRecord rec;
  BamRecordVector brv, brv2;
  size_t count = 0;
  while(br.GetNextRecord(rec) && count++ < 10000) 
    b.AddSequence(rec.Sequence().c_str(), rec.Qualities().c_str(), rec.Qname().c_str());
  
  b.Train();
  b.ErrorCorrect();

  float kcov = b.GetKCov();
  int   kmer = b.GetKMer();

  UnalignedSequenceVector v;

  v.clear();
  std::string seq, name;
  while (b.GetSequence(seq, name))
    v.push_back({name, seq});

  //std::ofstream filt("filt.fa");
  //for (auto& i : v) {
  //  filt << ">" << i.Name << std::endl << i.Seq << std::endl;
  //}
  //filt.close();

  FermiAssembler f;
  f.AddReads(v);
  f.DirectAssemble(kcov);

  // retrieve the contigs
  std::vector<std::string> contigs = f.GetContigs();

  std::ofstream cont("contigs.fa");
  size_t cc = 0;
  for (auto& i : f.GetContigs()) {
    ++cc;
    cont << ">" << cc << std::endl << i << std::endl;
  }
  cont.close();
  
}

BOOST_AUTO_TEST_CASE( header_check ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  // a BAM header check
  BOOST_CHECK_EQUAL(h.GetSequenceLength(0), 249250621);
  BOOST_CHECK_EQUAL(h.GetSequenceLength(3), 191154276);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("1"), 249250621);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("4"), 191154276);
  BOOST_CHECK_EQUAL(h.GetSequenceLength("d4"), -1);
  BOOST_CHECK_EQUAL(h.GetSequenceLength(10000), -1);
  
  BOOST_CHECK_EQUAL(h.GetHeaderSequenceVector().size(), h.NumSequences());
  BOOST_CHECK_EQUAL(h.GetHeaderSequenceVector().begin()->Length, 249250621);

  SeqLib::BamRecord rec;
  size_t count = 0;
  
  while(br.GetNextRecord(rec) && count++ < 10000) {
    
  }

}

BOOST_AUTO_TEST_CASE( merge ) {

  SeqLib::GRC grc;
  // add two more that we know of
  grc.add(SeqLib::GenomicRegion(23, 10,100));
  grc.add(SeqLib::GenomicRegion(23, 20,110));

  grc.add(SeqLib::GenomicRegion(2, 10,100));
  grc.add(SeqLib::GenomicRegion(2, 20,110));
  grc.add(SeqLib::GenomicRegion(2, 200,310));

  grc.MergeOverlappingIntervals();
  BOOST_CHECK_EQUAL(grc.size(), 3);
  BOOST_CHECK_EQUAL(grc[0].chr, 2);
  BOOST_CHECK_EQUAL(grc[1].chr, 2);
  BOOST_CHECK_EQUAL(grc[2].chr, 23);
  BOOST_CHECK_EQUAL(grc[2].pos2, 110);
  BOOST_CHECK_EQUAL(grc[2].pos1, 10);
}

BOOST_AUTO_TEST_CASE ( interval_queries ) {

  SeqLib::GRC grc;

  // create a large GRC
  for (int i = 0; i < 10; ++i) {
    int chr = rand() % 23;
    int pos = rand() % 10000;
    grc.add(SeqLib::GenomicRegion(chr, pos, pos + 100));
  }
  grc.MergeOverlappingIntervals();

  // add two more that we know of
  grc.add(SeqLib::GenomicRegion(23, 10,100));
  grc.add(SeqLib::GenomicRegion(23, 20,110));

  // create the interval tree
  grc.CreateTreeMap();

  SeqLib::GRC results = grc.FindOverlaps(SeqLib::GenomicRegion(23, 10, 100), true);

  for (auto& i : results)
    std::cerr << " GRC overlaps results " << i << std::endl;
  
  BOOST_CHECK_EQUAL(results.size(), 2);
  BOOST_CHECK_EQUAL(results[1].pos2, 100);

  grc.MergeOverlappingIntervals();
  grc.CreateTreeMap();

  for(auto& r : grc)
    std::cerr << r << std::endl;

  std::vector<int32_t> q, s;
  results = grc.FindOverlaps(grc, q, s, true);

  std::cerr << " results.size " << results.size() << " Input size " << grc.size() << std::endl;
  BOOST_CHECK_EQUAL(results.size(), grc.size());
  BOOST_CHECK_EQUAL(results.TotalWidth(), grc.TotalWidth());
  
}

BOOST_AUTO_TEST_CASE( json_parse_from_file ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  
  std::string rules = "{\"global\" : {\"!anyflag\" : 1536}, \"\" : { \"rules\" : [{\"ic\" : true}, {\"clip\" : 5}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  ReadFilterCollection rfc(rules, br.Header());

  std::cerr << rfc << std::endl;

  SeqLib::BamRecord rec;
  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) {
    if (!rfc.isValid(rec))
      continue;
    // test global flag rule
    if ( (rec.QCFailFlag() || rec.DuplicateFlag())) {
      std::cerr << rec << std::endl;
      assert(false);
    }
  }

  /// direct from string
  ReadFilterCollection rfc2(rules, br.Header());

  while(br.GetNextRecord(rec) && count++ < 10000) {
    if (!rfc.isValid(rec))
      continue;
    // test global flag rule
    if ( (rec.QCFailFlag() || rec.DuplicateFlag())) {
      std::cerr << rec << std::endl;
      assert(false);
    }
  }

  // check that a bad key throws error
  //rules = "{\"global\" : {\"!anyflagf\" : 1536}, \"\" : { \"rules\" : [{\"ic\" : true}, {\"clip\" : 5}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  //BOOST_CHECK_THROW(ReadFilterCollection rfc2(rules, br.Header()), std::invalid_argument);
  // bad value, expected int
  //rules = "{\"global\" : {\"!anyflag\" : \"BAD\"}, \"\" : { \"rules\" : [{\"ic\" : true}, {\"clip\" : 5}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  //BOOST_CHECK_THROW(ReadFilterCollection rfc3(rules, br.Header()), std::invalid_argument);
  // bad JSON itself
  rules = "{\"global\" : \"!anyflag\" : 1536}, \"\" : { \"rules\" : [{\"ic\" : true}, {\"clip\" : 5}, {\"ins\" : true}, {\"del\" : true}, {\"mapped\": true , \"mate_mapped\" : false}, {\"mate_mapped\" : true, \"mapped\" : false}]}}";  
  BOOST_CHECK_THROW(ReadFilterCollection rfc4(rules, br.Header()), std::invalid_argument);
  // bad value, expected range
  rules = "{\"global\" : {\"!anyflag\" : 1536}, \"\" : { \"rules\" : [{\"isize\" : \"BAD\"}]}}";  
  BOOST_CHECK_THROW(ReadFilterCollection rfc4(rules, br.Header()), std::invalid_argument);

    
}

BOOST_AUTO_TEST_CASE( sw_alignment ) {

  const std::string ref = "ACTGCGAGCGACTAGCTCGTAGCTAGCTAGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCGATCGCTACGAC";
  const std::string seq = "ACTGCGAGCGACTAGCTCGTAGCTAGCTAGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCGATCGCTACGAC";
  //const std::string seq =                "CTCGTAGCTAGCTGCTAGCTAGTGACTGCGGGCGATCATCGATCTTTTATTATCGCG";
  const SeqLib::GenomicRegion gr(0,0,0);
  SeqLib::BamRecord b("test_name", seq, ref, &gr);
  
  std::cerr << " SMITH WATERMAN " << std::endl;
  std::cerr << b << std::endl;
}

BOOST_AUTO_TEST_CASE( read_filter_1 ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::GRC g;
  g.add(SeqLib::GenomicRegion(h.Name2ID("X"), 1100000, 1800000));
  g.CreateTreeMap();

  // make a new rule set
  ReadFilterCollection rfc;

  // make a new filter region
  ReadFilter rf;

  // add an isize rule on whole-genome
  AbstractRule ar;
  ar.isize = Range(200, 600, false); // 200 to 600, not inverted
  ar.mapq  = Range(10, 50, false); // 200 to 600, not inverted
  ar.nm    = Range(1, 1, false); // 200 to 600, not inverted
  rf.AddRule(ar);

  rf.setRegions(g);

  // add to the filter collection
  rfc.AddReadFilter(rf);

  SeqLib::GRC gback = rfc.getAllRegions();
  BOOST_CHECK_EQUAL(gback.size(), g.size());
  for (size_t i = 0; i < gback.size(); ++i)
    assert(g[i] == gback[i]);
  

  // display
  std::cerr << br.PrintRegions() << std::endl;

  // read / filter the reads
  SeqLib::BamRecord rec;
  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) {

    if (!rfc.isValid(rec))
      continue;

    // test isize rule
    if (!(rec.FullInsertSize() >= 200 || rec.FullInsertSize() <= 600)) {
      std::cerr << rec.FullInsertSize() << std::endl;
      assert(false);
    }
    // test mapq rule
    if (!(rec.MapQuality() >= 10 || rec.MapQuality() <= 50)) {
      std::cerr << rec.MapQuality() << std::endl;
      assert(false);
    }
    // test nm rule
    int32_t nm;
    rec.GetIntTag("NM", nm);
    if (nm == 1) 
      assert(false);
  }
}

BOOST_AUTO_TEST_CASE ( fermi_add_reads ) {

  SeqLib::FermiAssembler f;
  
  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamRecord r;
  size_t count = 0;
  while (br.GetNextRecord(r) && count++ < 1000) 
    f.AddRead(r);
  while (br.GetNextRecord(r) && count++ < 2000)   
    f.AddRead(SeqLib::UnalignedSequence(r.Qname(), r.Sequence(), r.Qualities()));
  
  f.CorrectReads();
  f.PerformAssembly();
  std::vector<std::string> out = f.GetContigs();
  
}

BOOST_AUTO_TEST_CASE ( seq_utils ) {
  
  // add commas
  BOOST_CHECK_EQUAL(SeqLib::AddCommas(1),"1");
  BOOST_CHECK_EQUAL(SeqLib::AddCommas(1000000),"1,000,000");
  
  // percent calc
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(10,100), 10);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(7,8), 87);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(9,10), 90);
  BOOST_CHECK_EQUAL(SeqLib::percentCalc(2,3), 66);

  // scrub string
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "chr"), "1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", ""), "chr1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "dd"), "chr1");
  BOOST_CHECK_EQUAL(SeqLib::scrubString("chr1", "1"), "chr");

}

BOOST_AUTO_TEST_CASE( bam_record ) {

  // get a record
  SeqLib::BamReader br, br2;

  // try multiple open
  std::vector<std::string> bs = {"test_data/small.bam", "test_data/small.bam"};
  BOOST_CHECK(!br2.Open(bs)); // should be false, no dups
  bs = {"test_data/small.bam", "test_data/small.cram"};
  BOOST_CHECK(br.Open(bs)); // should be true
  SeqLib::BamRecord r;
  
  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  br.GetNextRecord(r);
  
  BOOST_CHECK_EQUAL(r.AsGenomicRegion().chr, 0);
  BOOST_CHECK_EQUAL(r.AsGenomicRegion().pos1,9995);
  BOOST_CHECK_EQUAL(r.AsGenomicRegion().pos2,10075);
  BOOST_CHECK_EQUAL(r.AsGenomicRegion().strand,'-');

  BOOST_CHECK_EQUAL(r.AsGenomicRegionMate().chr, 15);
  BOOST_CHECK_EQUAL(r.AsGenomicRegionMate().pos1,67300983);
  BOOST_CHECK_EQUAL(r.AsGenomicRegionMate().pos2,67301134);
  BOOST_CHECK_EQUAL(r.AsGenomicRegionMate().strand,'-');

  BOOST_CHECK_EQUAL(std::floor(r.MeanPhred()), 15);

  BOOST_CHECK_EQUAL(r.CountNBases(), 0);

  r.SetQname("testq");
  BOOST_CHECK_EQUAL(r.Qname(), "testq");

  const std::string s = "ACTGCTAGCTAGCTACTCTGCTACTATATTAGCGCGCATTCGC";
  r.SetSequence(s);
  BOOST_CHECK_EQUAL(r.Sequence(), s);
  
  r.SmartAddTag("ST", "1");
  r.SmartAddTag("ST", "3");
  r.SmartAddTag("ST", "5");
  r.SmartAddTag("S2", "5");
    
  BOOST_CHECK_EQUAL(r.GetSmartIntTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartIntTag("ST").at(2), 5);
  BOOST_CHECK_EQUAL(r.GetSmartDoubleTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartDoubleTag("ST").at(2), 5.0);
    
  BOOST_CHECK_EQUAL(r.GetSmartStringTag("ST").size(), 3);
  BOOST_CHECK_EQUAL(r.GetSmartStringTag("ST")[1], "3");

  BOOST_CHECK_EQUAL(r.GetSmartDoubleTag("S2").at(0), 5.0);
  BOOST_CHECK_EQUAL(r.GetSmartIntTag("S2").at(0), 5);
}

BOOST_AUTO_TEST_CASE( fermi_assemble ) {

  SeqLib::FermiAssembler f;

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamRecord r;

  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  while(br.GetNextRecord(r) && count++ < 1000) {
    brv.push_back(r);
  }
  
  f.AddReads(brv);

  f.CorrectReads();

  SeqLib::UnalignedSequenceVector reads = f.GetSequences();
  BOOST_CHECK_EQUAL(reads.size(), brv.size());

  for (int i = 0; i < reads.size(); ++i) {
    if (brv[i].Sequence() != reads[i].Seq) {
      std::cerr << "************" << std::endl;
      std::cerr << brv[i].Sequence() << std::endl;
      std::cerr << reads[i].Seq << std::endl;
    }
  }
    
  // peform the assembly
  std::cerr << "...performing assembly" << std::endl;
  f.PerformAssembly();

  // retrieve the contigs
  std::vector<std::string> contigs = f.GetContigs();

}


BOOST_AUTO_TEST_CASE( bam_header_stdout ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  std::cout << h.AsString() << std::endl;
}

BOOST_AUTO_TEST_CASE( bam_header_name2id ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  BOOST_CHECK_EQUAL(h.Name2ID("2"), 1);  
  BOOST_CHECK_EQUAL(h.Name2ID("23"), -1);  

}

BOOST_AUTO_TEST_CASE( bam_header_id2name ) {
  
  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();
  
  BOOST_CHECK_EQUAL(h.IDtoName(2), "3");
  BOOST_CHECK_THROW(h.IDtoName(100), std::out_of_range);
  BOOST_CHECK_THROW(h.IDtoName(-1), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::BamHeader().IDtoName(1), std::out_of_range);
}

BOOST_AUTO_TEST_CASE( genomic_ranges_string_constructor) {
  
  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();
  
  const std::string in = "2:1,000,000-2,000,000";
  SeqLib::GenomicRegion gr(in, h);
  BOOST_CHECK_EQUAL(gr.chr, 1);
  BOOST_CHECK_EQUAL(gr.pos1, 1000000);
  BOOST_CHECK_EQUAL(gr.pos2, 2000000);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(in, SeqLib::BamHeader()), std::invalid_argument);

  BOOST_CHECK_EQUAL(gr.ChrName(h), "2");
  BOOST_CHECK_EQUAL(gr.ChrName(SeqLib::BamHeader()), "2");
  gr.chr = 1000;
  BOOST_CHECK_THROW(gr.ChrName(h), std::invalid_argument);

}

BOOST_AUTO_TEST_CASE( genomic_region_less_than ) {

  SeqLib::GenomicRegion gr1(0, 1, 2);
  SeqLib::GenomicRegion gr2(1, 1, 2);
  SeqLib::GenomicRegion gr3(1, 2, 2);
  SeqLib::GenomicRegion gr4(1, 6, 6);

  BOOST_CHECK(gr1 < gr2);
  BOOST_CHECK(gr2 > gr1);
  BOOST_CHECK(!(gr1 > gr2));

  BOOST_CHECK(gr2 < gr3);
  BOOST_CHECK(gr3 > gr2);
  BOOST_CHECK(!(gr2 > gr3));

  BOOST_CHECK(gr3 < gr4);
  BOOST_CHECK(!(gr4 == gr3));
  BOOST_CHECK(!(gr3 > gr4));
  BOOST_CHECK(gr4 > gr3);

  BOOST_CHECK(!(gr1 < gr1));
  BOOST_CHECK(!(gr1 > gr1));

  BOOST_CHECK(!(gr1 != gr1));
  BOOST_CHECK(gr2 != gr1);
  BOOST_CHECK(gr3 != gr1);
  BOOST_CHECK(gr4 != gr3);

  BOOST_CHECK(gr1 >= gr1);
  BOOST_CHECK(gr2 >= gr2);
  BOOST_CHECK(gr3 >= gr3);
  BOOST_CHECK(gr4 >= gr4);

  BOOST_CHECK(gr1 <= gr1);
  BOOST_CHECK(gr2 <= gr2);
  BOOST_CHECK(gr3 <= gr3);
  BOOST_CHECK(gr4 <= gr4);

  BOOST_CHECK(gr1 <= gr2);
  BOOST_CHECK(gr2 >= gr1);

  BOOST_CHECK(gr2 <= gr3);
  BOOST_CHECK(gr3 >= gr2);

}

BOOST_AUTO_TEST_CASE( genomic_region_distance ) {

  SeqLib::GenomicRegion gr1(0, 10, 100);
  SeqLib::GenomicRegion gr2(0, 10, 200);
  SeqLib::GenomicRegion gr3(1, 10, 100);
  SeqLib::GenomicRegion gr4(0, 100, 100);

  BOOST_CHECK_EQUAL(gr1.DistanceBetweenEnds(gr3), -1);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenEnds(gr1), 0);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenEnds(gr2), 100);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenEnds(gr4), 0);

  BOOST_CHECK_EQUAL(gr1.DistanceBetweenStarts(gr3), -1);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenStarts(gr1), 0);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenStarts(gr2), 0);
  BOOST_CHECK_EQUAL(gr1.DistanceBetweenStarts(gr4), 90);

}

BOOST_AUTO_TEST_CASE( small_trie_from_file) {

  AbstractRule ar;
  const bool inverted = false;
  ar.addMotifRule("test_data/motif.txt", inverted);

  ReadFilterCollection rfc;
  ReadFilter rf;
  rf.AddRule(ar);
  rfc.AddReadFilter(rf);

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");

  SeqLib::BamRecord rec;
  bool rule;
  size_t count = 0;
  while (br.GetNextRecord(rec) && count++ < 1000){
  }
  
}

BOOST_AUTO_TEST_CASE( large_trie ) {

  const std::string dictionary = "ACTG";
  
  const int string_size = 20;
  const int string_count = 10000;

  SeqLib::Filter::AhoCorasick aho;

  std::vector<std::string> k;

  std::cerr << "...generating key" << std::endl;
  for (int i = 0; i < string_count; ++i) {
    char* c = (char*) malloc(string_size + 1);
    for (int j = 0; j < string_size; ++j)
      c[j] = dictionary.at(rand() % 4);
    c[string_size] = '\0';
    k.push_back(std::string(c));
    free(c);
  }
  std::cerr << "...done with key" << std::endl;

  std::cerr << "...generating trie" << std::endl;
  for (auto& i : k)
    aho.AddMotif(i);
  std::cerr << "...done generating trie" << std::endl;

  std::cerr << "...querying trie" << std::endl;
  auto result = aho.aho_trie->parse_text(k[0]);
  std::cerr << "...querying trie fast" << std::endl;  
  for (int i = 0; i < string_count; ++i) {
    //if (i % 20000 == 0)
    //  std::cerr << "... " << i << std::endl;
    auto result = aho.aho_trie->parse_text(k[i]);
  }
    
}

BOOST_AUTO_TEST_CASE( genomic_region_constructors ) {

  // GenomicRegion Constructors
  SeqLib::GenomicRegion gr(0, 0, 10, '+');
  BOOST_CHECK_EQUAL(gr.Width(), 11);

  SeqLib::GenomicRegion gr_empty;
  BOOST_TEST(gr_empty.IsEmpty());

  SeqLib::GenomicRegion gr2("chrX", "0", "10", SeqLib::BamHeader());
  BOOST_CHECK_EQUAL(gr2.Width(), 11);
  BOOST_CHECK_EQUAL(gr2.chr, 22);

  SeqLib::GenomicRegion gr3("X", "0", "10", SeqLib::BamHeader());
  BOOST_TEST(gr2 == gr3);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion gr3("X", "a", "10", SeqLib::BamHeader()), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::GenomicRegion gr3("X", "1000000000000000000000000000000000000000000000000000000000000000000000000000000", "10", SeqLib::BamHeader()), std::out_of_range);

  BOOST_CHECK_EQUAL(gr.DistanceBetweenStarts(gr2), -1);
  BOOST_CHECK_EQUAL(gr2.DistanceBetweenStarts(gr), -1);

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  BOOST_CHECK_EQUAL(SeqLib::GenomicRegion("X","1","100", br.Header()).chr, 22);

  // check negative inputs
  SeqLib::GenomicRegion grn(-1,-11,-10);
  BOOST_CHECK_EQUAL(grn.chr, -1);
  BOOST_CHECK_EQUAL(grn.pos1, -11);
  BOOST_CHECK_EQUAL(grn.pos2, -10);

  // check strand constructions
  SeqLib::GenomicRegion gra(0,0,0);
  SeqLib::GenomicRegion grb(0,10000,10001, '+');
  SeqLib::GenomicRegion grc(0,0,3, '-');
  BOOST_CHECK_EQUAL(gra.strand, '*');
  BOOST_CHECK_EQUAL(grb.strand, '+');
  BOOST_CHECK_EQUAL(grc.strand, '-');

  // check point string
  BOOST_CHECK_EQUAL(grb.PointString(), "1:10,000(+)");

  // check pretty string
  std::stringstream ss;
  ss << grb;
  BOOST_CHECK_EQUAL(ss.str(), "1:10,000-10,001(+)");

}

BOOST_AUTO_TEST_CASE( genomic_region_bad_inputs ) {

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(0, 10, 9), std::invalid_argument);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion(0,0,0,'P'), std::invalid_argument);

}

// BOOST_AUTO_TEST_CASE( genomic_region_random ) {

//   SeqLib::GenomicRegion gr; 
//   std::srand(42);
//   gr.random();
//   BOOST_CHECK_EQUAL(gr.pointString(), "9:69,477,830(*)");
  
// }

BOOST_AUTO_TEST_CASE( genomic_region_range_operations ) {

  SeqLib::GenomicRegion gr(0,1,10);
  SeqLib::GenomicRegion gr2(0,1,11);
  gr.Pad(3);
  gr2.Pad(-3);
  BOOST_CHECK_EQUAL(gr.pos1,-2);
  BOOST_CHECK_EQUAL(gr.pos2,13);
  BOOST_CHECK_EQUAL(gr2.pos1,4);
  BOOST_CHECK_EQUAL(gr2.pos2,8);

  BOOST_CHECK_THROW(gr.Pad(-10), std::out_of_range);

}

BOOST_AUTO_TEST_CASE( genomic_check_overlaps ) {

  SeqLib::GenomicRegion gr1(0, 0, 10, '+');
  SeqLib::GenomicRegion gr2(1, 0, 10, '+');

  SeqLib::GenomicRegion gr3(0, 10, 20, '+');
  SeqLib::GenomicRegion gr4(1, 4, 10, '+');

  SeqLib::GenomicRegion gr5(1, 11, 12, '+');

  // partial overlaps should be one
  BOOST_CHECK_EQUAL(gr1.GetOverlap(gr3), 1);

  // argument contained gets 2
  BOOST_CHECK_EQUAL(gr2.GetOverlap(gr4), 2);

  // object contained gets 3 
  BOOST_CHECK_EQUAL(gr4.GetOverlap(gr2), 3);

  // same chr, no overlap
  BOOST_CHECK_EQUAL(gr4.GetOverlap(gr5), 0);
  BOOST_CHECK_EQUAL(gr5.GetOverlap(gr4), 0);

}

BOOST_AUTO_TEST_CASE( bwa_wrapper ) {

  SeqLib::BWAWrapper bwa;

  // Set some options
  bwa.SetGapOpen(32);
  bwa.SetGapExtension(1);
  bwa.SetMismatchPenalty(18);
  bwa.SetAScore(2);
  bwa.SetZDropoff(100);
  bwa.Set3primeClippingPenalty(5);
  bwa.Set5primeClippingPenalty(5);
  bwa.SetBandwidth(1000);
  bwa.SetReseedTrigger(1.5);

  BOOST_CHECK_THROW(bwa.SetGapOpen(-1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetGapExtension(-1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetMismatchPenalty(-18), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetAScore(-2), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetZDropoff(-100), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.Set3primeClippingPenalty(-5), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.Set5primeClippingPenalty(-5), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetBandwidth(-1000), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.SetReseedTrigger(-1.5), std::invalid_argument);

  // no index loaded yet
  BOOST_CHECK_THROW(bwa.ChrIDToName(1), std::runtime_error);

  // load a test index
  BOOST_TEST(SeqLib::read_access_test(TREF));
  bwa.LoadIndex(TREF);

  BOOST_CHECK_EQUAL(bwa.NumSequences(), 2);

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref1");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref2");
  BOOST_CHECK_THROW(bwa.ChrIDToName(2), std::out_of_range);

  BOOST_CHECK(!bwa.LoadIndex("test_data/small.bam"));

  SeqLib::BamHeader hh = bwa.HeaderFromIndex();
  BOOST_CHECK_EQUAL(hh.NumSequences(), 2);

  // error check the index construction
  SeqLib::UnalignedSequenceVector usv_bad1, usv_bad2, usv;;
  usv_bad1.push_back(SeqLib::UnalignedSequence("ref1","ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC", std::string()));
  usv_bad1.push_back(SeqLib::UnalignedSequence("ref4", std::string(), std::string()));
  usv_bad1.push_back(SeqLib::UnalignedSequence("ref5","CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT", std::string()));
  usv_bad2.push_back(SeqLib::UnalignedSequence(std::string(), "ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC", std::string()));
  usv_bad2.push_back(SeqLib::UnalignedSequence("ref4","ACCATCGCAGCAGCTATCTATTATATCGGCAGCATCTAGC", std::string()));
  usv_bad2.push_back(SeqLib::UnalignedSequence("ref5","CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT", std::string()));
  BOOST_CHECK_THROW(bwa.ConstructIndex(usv_bad1), std::invalid_argument);
  BOOST_CHECK_THROW(bwa.ConstructIndex(usv_bad2), std::invalid_argument);
  
  // construct a normal index
  usv.push_back(SeqLib::UnalignedSequence("ref3","ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCGATCGATCGATCGTAGC", std::string()));
  usv.push_back(SeqLib::UnalignedSequence("ref4","CTACTTTATCATCTACACACTGCCTGACTGCGGCGACGAGCGAGCAGCTACTATCGACT", std::string()));
  usv.push_back(SeqLib::UnalignedSequence("ref5","CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT", std::string()));
  usv.push_back(SeqLib::UnalignedSequence("ref6","TATCTACTGCGCGCGATCATCTAGCGCAGGACGAGCATC" + std::string(100,'N') + "CGATCGTTATTATCGAGCGACGATCTACTACGT", std::string()));

  bwa.ConstructIndex(usv);

  BOOST_CHECK_EQUAL(bwa.NumSequences(), 4);
  bwa.ChrIDToName(1);

  BOOST_CHECK_THROW(bwa.ChrIDToName(-1), std::out_of_range);
  BOOST_CHECK_THROW(bwa.ChrIDToName(10000), std::out_of_range);

  std::cout << bwa.ChrIDToName(2) << std::endl;

  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref3");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref4");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(2), "ref5");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(3), "ref6");
  BOOST_CHECK_THROW(bwa.ChrIDToName(4), std::out_of_range);

  // write the index
  BOOST_CHECK(bwa.WriteIndex(OREF));

  // write the fasta
  std::ofstream os;
  os.open(OREF);
  os          << ">" << usv[0].Name << std::endl << usv[0].Seq <<
    std::endl << ">" << usv[1].Name << std::endl << usv[1].Seq << 
    std::endl << ">" << usv[2].Name << std::endl << usv[2].Seq << 
    std::endl << ">" << usv[3].Name << std::endl << usv[3].Seq << 
    std::endl;
  os.close();

  // read it back
  BOOST_CHECK(bwa.LoadIndex(OREF));

  // check that its good
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(0), "ref3");
  BOOST_CHECK_EQUAL(bwa.ChrIDToName(1), "ref4");
  
  // try aligning a sequence
  std::cerr << "...aligning sequences" << std::endl;
  SeqLib::BamRecordVector brv, brv2;
  bool hardclip = false;
  bwa.AlignSequence("ACATGGCGAGCACTTCTAGCATCAGCTAGCTACGATCG", "name", brv, 0.9, hardclip, 1);
  // reverse complement
  bwa.AlignSequence("CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGC", "name", brv2, 0.9, hardclip, 2);

  BOOST_CHECK_EQUAL(brv[0].Qname(), "name");
  BOOST_CHECK_EQUAL(brv[0].ChrID(), 2);
  BOOST_CHECK_EQUAL(brv[0].Sequence(), "CGATCGTAGCTAGCTGATGCTAGAAGTGCTCGCCATGT");
  std::cerr << " brv[0].GetCigar() " << brv[0].GetCigar() << std::endl;
  BOOST_CHECK_EQUAL(brv[0].GetCigar()[0].Type(), 'M');
  BOOST_CHECK_EQUAL(brv[0].GetCigar()[0].Length(), 38);

  // check from iterator
  SeqLib::Cigar ccc = brv[0].GetCigar();
  
  //SeqLib::Cigar::const_iterator f = brv[0].GetCigar().begin();
  BOOST_CHECK_EQUAL(ccc.begin()->Length(), 38);

  // check that it got both alignments
  BOOST_CHECK_EQUAL(brv2.size(), 2);

  // print info 
  std::cerr << bwa << std::endl;
}

BOOST_AUTO_TEST_CASE( bam_reader ) {

  // print empty read
  std::cerr << SeqLib::BamRecord() << std::endl;

  SeqLib::BamReader bw;
  bw.Open(SBAM);

  // open index
  bw.SetRegion(SeqLib::GenomicRegion(22, 1000000, 1001000));

  // make a set of locations
  SeqLib::GRC grc;
  for (size_t i = 0; i < 24; ++i)
    grc.add(SeqLib::GenomicRegion(i, 1, 100));

  // set regions
  bw.SetMultipleRegions(grc);

  // write index of new bam
  // should print a warning since no write bam is specified
  //bw.BuildIndex();

  // open an output BAM
  //bw.OpenWriteBam(OBAM);

  // set tags to strip
  //bw.setStripTags("OQ,BI");

  // loop through and grab some reads
  SeqLib::BamRecord r;
  size_t count = 0;
  while (bw.GetNextRecord(r)) {
    //if (++count % 10 == 0)
    //  bw.WriteRecord(r);
  }
  
  // display info about BAM
  std::cerr << bw << std::endl;

  // write index of new bam
  //bw.BuildIndex();

  // reset the walker
  bw.Reset();

  // set a smaller region
  bw.SetRegion(SeqLib::GenomicRegion(0, 1, 100));  
  std::cerr << bw << std::endl;

  // write as a cram
  //bw.OpenWriteBam(OCRAM);
    
  //
  //bw.setCram(OCRAM, HGREF);

  // print cram writer
  //std::cerr << bw << std::endl;
  // write the CRAM
  //while (bw.GetNextRecord(r, rule)) {
  //  if (++count % 10 == 0) {
  //    std::cerr << count << std::endl;
  //    bw.WriteRecord(r);
  //  }
  //}

}

BOOST_AUTO_TEST_CASE( set_qualities ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");


  SeqLib::BamRecord r;
  while (br.GetNextRecord(r)) {
    r.SetQualities("", 0);
    BOOST_CHECK_EQUAL(r.Qualities().at(0), '!');
    r.SetQualities(std::string(r.Length(), '#'), 33);
    BOOST_CHECK_EQUAL(r.Qualities(), std::string(r.Length(), '#'));    
    BOOST_CHECK_THROW(r.SetQualities(std::string(8, '#'), 0), std::invalid_argument);
    break;
  }
}

BOOST_AUTO_TEST_CASE( header_constructor ) {

  SeqLib::HeaderSequenceVector hsv;
  hsv.push_back(SeqLib::HeaderSequence("1", 1000));
  hsv.push_back(SeqLib::HeaderSequence("chr2", 1200));
  SeqLib::BamHeader hdr(hsv);
  
}

BOOST_AUTO_TEST_CASE( overlapping_coverage ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamRecordVector brv;
  size_t count = 0;
  SeqLib::BamRecord r;
  while(br.GetNextRecord(r) && ++count < 4) {
    std::cout << " r " << r << std::endl;
    brv.push_back(r);
  }
  BOOST_CHECK_EQUAL(brv[0].OverlappingCoverage(brv[2]), 78);

}

BOOST_AUTO_TEST_CASE( gr_chr_region_set) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  

  SeqLib::GenomicRegion gr("1", br.Header());
  BOOST_CHECK_EQUAL(gr.chr, 0);
  BOOST_CHECK_EQUAL(gr.pos2, 249250621);
  BOOST_CHECK_EQUAL(gr.pos1, 1);

  BOOST_CHECK_THROW(SeqLib::GenomicRegion gr2("-1", br.Header()), std::invalid_argument);

}

BOOST_AUTO_TEST_CASE( sequtils ) {

  std::string seq = "actgACGTnTCN";

  SeqLib::rcomplement(seq);
  
  BOOST_CHECK_EQUAL(seq, "NGAnACGTcagt");

}

BOOST_AUTO_TEST_CASE( bam_write ) {


  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::BamRecord rec;

  // empty constructor
  SeqLib::BamWriter w;
  
  BOOST_CHECK(!w.WriteHeader());
  BOOST_CHECK(!w.Close());
  BOOST_CHECK(!w.BuildIndex());
  BOOST_CHECK(!w.WriteRecord(rec));

  w.Open("tmp_out.bam");

  // check that set CRAM fails
  BOOST_CHECK(!w.SetCramReference("dummy")); 
  BOOST_CHECK(!w.WriteHeader());

  w.SetHeader(h);

  w.WriteHeader();

  size_t count = 0;

  while(br.GetNextRecord(rec) && count++ < 10000) 
    w.WriteRecord(rec);

  BOOST_CHECK(!w.BuildIndex());
  w.Close();

  w.BuildIndex();

  // check that write header now fails
  BOOST_CHECK(!w.WriteHeader());

  // check that set CRAM fails
  BOOST_CHECK(!w.SetCramReference("badref"));

  // print some info
  std::cerr << w << std::endl;

}

BOOST_AUTO_TEST_CASE( bam_record_more ) {

  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamHeader h = br.Header();

  SeqLib::BamRecord rec;
  size_t count = 0;
  
  while(br.GetNextRecord(rec) && count++ < 100) {
    rec.ClearSeqQualAndTags();
    assert(rec.Sequence().empty());
    assert(rec.Qualities().empty());
    //int32_t nm;
    //rec.GetIntTag("NM", nm);
    //assert(!nm);
    std::string xa;
    rec.GetZTag("XA", xa);
    BOOST_CHECK_EQUAL(xa.empty(), rec.CountBWASecondaryAlignments()==0);
    rec.CountBWAChimericAlignments();
  }

  br.Reset();

  SeqLib::Filter::ReadFilterCollection rf;

}

BOOST_AUTO_TEST_CASE( bam_record_manipulation ) {

  SeqLib::Cigar cig;

  // manually construct a cigar
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('I', 1));
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('D', 1));
  cig.add(SeqLib::CigarField('M', 10));
  cig.add(SeqLib::CigarField('S', 10));

  // check that coversion to cigar data strutur (uint32_t) worked
  SeqLib::CigarField cm('M', 1);
  SeqLib::CigarField ci('I', 1);
  SeqLib::CigarField cd('D', 1);
  SeqLib::CigarField cn('N', 1);
  SeqLib::CigarField cs('S', 1);
  SeqLib::CigarField ch('H', 1);
  SeqLib::CigarField cp('P', 1);
  SeqLib::CigarField ce('=', 1);
  SeqLib::CigarField cx('X', 1);
  SeqLib::CigarField cb('B', 1);

  BOOST_CHECK_EQUAL(cm.Type(), 'M');
  BOOST_CHECK_EQUAL(ci.Type(), 'I');
  BOOST_CHECK_EQUAL(cd.Type(), 'D');
  BOOST_CHECK_EQUAL(cn.Type(), 'N');
  BOOST_CHECK_EQUAL(cs.Type(), 'S');
  BOOST_CHECK_EQUAL(ch.Type(), 'H');
  BOOST_CHECK_EQUAL(cp.Type(), 'P');
  BOOST_CHECK_EQUAL(ce.Type(), '=');
  BOOST_CHECK_EQUAL(cx.Type(), 'X');

  // check invalid constructions
  BOOST_CHECK_THROW(SeqLib::CigarField('L', 1), std::invalid_argument);

  // make a sequence
  const std::string seq = std::string(10, 'A') + std::string(1, 'T') + std::string(10, 'C') + std::string(10, 'G') + std::string(10, 'A');

  // check   
  BOOST_CHECK_EQUAL(cig.NumQueryConsumed(), 41);
  BOOST_CHECK_EQUAL(cig.NumReferenceConsumed(), 31);

  std::stringstream ss;
  ss << cig;

  // cigar from string
  SeqLib::Cigar cig2(ss.str());

  // check that the string from / to are consistent
  assert(cig == cig2);
  assert(!(cig != cig2));
  for (int i = 0; i < cig.size(); ++i)
    assert(cig[i] == cig2[i]);
  for (int i = 0; i < cig.size(); ++i)
    assert(!(cig[i] != cig2[i]));

  // manually make a read
  SeqLib::GenomicRegion gr_wrong(0, 100, 131); 
  SeqLib::GenomicRegion gr(0, 100, 130); 
  
  BOOST_CHECK_THROW(SeqLib::BamRecord("dumname", seq, &gr_wrong, cig), std::invalid_argument);
  BOOST_CHECK_THROW(SeqLib::BamRecord("dumname", seq + "A", &gr, cig), std::invalid_argument);

  SeqLib::BamRecord br("dumname", seq, &gr, cig);

  BOOST_CHECK_EQUAL(br.Sequence(), seq);
  BOOST_CHECK_EQUAL(br.GetCigar(), cig);
  BOOST_CHECK_EQUAL(br.Qname(), "dumname");
  BOOST_CHECK_EQUAL(br.Position(), 100);
  BOOST_CHECK_EQUAL(br.Length(), 41);
  BOOST_CHECK_EQUAL(br.ChrID(), 0);

}

BOOST_AUTO_TEST_CASE( change_bam_record ) {

  // get a record
  SeqLib::BamReader br;
  br.Open("test_data/small.bam");
  SeqLib::BamRecord r;
  
  SeqLib::BamRecordVector brv;
  
  size_t count = 0;
  br.GetNextRecord(r);

  SeqLib::Cigar c = r.GetCigar();
  std::cerr << c << std::endl;

  // try replace with cigar of same size
  SeqLib::Cigar c2;
  c2.add(SeqLib::CigarField('S', 101));
  r.SetCigar(c2);
  std::cerr << r << std::endl;

  // try replace with new cigar
  SeqLib::Cigar c3;
  c3.add(SeqLib::CigarField('S', 10));
  c3.add(SeqLib::CigarField('M', 91));
  r.SetCigar(c3);
  std::cerr << r << std::endl;

  const std::string new_seq = "ACTGGACTACAC";

  r.SetSequence(new_seq);
  std::cerr << r << std::endl;

  r.SetQname("dummy_qname");
  std::cerr << r << std::endl;

}

BOOST_AUTO_TEST_CASE( stdinput ) {

#ifdef RUN_STDIN
  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("-"); 

  // write it back out
  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << " STDIN " << r << std::endl;
  }
#endif
}

BOOST_AUTO_TEST_CASE( cramin ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.cram"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << "CRAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( cramin_new_ref ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.cram");

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 10) {
    std::cerr << "CRAM " << r << std::endl;
  }

  b.Reset();
  
  // should fail
  

}


BOOST_AUTO_TEST_CASE( bamin ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.bam"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << "BAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( samin ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.sam"); 

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    std::cerr << "SAM " << r << std::endl;
  }
}

BOOST_AUTO_TEST_CASE( bamout ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::BAM);
  //SeqLib::BamWriter w;
  w.Open("tmp_out.bam");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  
}

BOOST_AUTO_TEST_CASE( samout ) {

  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("tmp_out.sam");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  b.Close();
}


BOOST_AUTO_TEST_CASE( cramout ) {

  SeqLib::BamReader b;
  b.Open("test_data/small.cram"); 

  SeqLib::BamWriter w(SeqLib::CRAM);
  w.Open("tmp_out.cram");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
  
}

BOOST_AUTO_TEST_CASE( samout_to_stdout ) {

#ifdef RUN_SAM_STDOUT
  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::SAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
#endif
}

BOOST_AUTO_TEST_CASE( bamout_to_stdout ) {

  //
  // dont actually run every time
  // too much stdout-ing
  //

#ifdef RUN_BAM_STDOUT
  // read a BAM from stdin
  SeqLib::BamReader b;
  b.Open("test_data/small.sam"); 

  SeqLib::BamWriter w(SeqLib::BAM);
  w.Open("-");
  w.SetHeader(b.Header());
  w.WriteHeader();

  SeqLib::BamRecord r;
  size_t count = 0;
  while(b.GetNextRecord(r) && count++ < 1) {
    w.WriteRecord(r);
  }
  w.Close();
#endif
  
}

BOOST_AUTO_TEST_CASE( bam_poly ) {

  SeqLib::BamReader r;
  
  BOOST_CHECK(r.Open("test_data/small.bam"));
  BOOST_CHECK(r.Open("test_data/small.cram"));

  BOOST_CHECK(r.SetRegion(SeqLib::GenomicRegion(r.Header().Name2ID("X"),1001000, 1001100)));
  BOOST_CHECK(!r.SetRegion(SeqLib::GenomicRegion(1000, 1001000, 1001100))); // should fail

  SeqLib::BamWriter w(SeqLib::BAM);
  w.Open("tmp_out_poly.bam");
  w.SetHeader(r.Header());
  w.WriteHeader();

  SeqLib::BamRecord rec;
  while(r.GetNextRecord(rec)) {
    w.WriteRecord(rec);
  }

  BOOST_CHECK(r.Reset("test_data/small.bam"));
  BOOST_CHECK(!r.Reset("dum"));

  BOOST_CHECK(r.Close("test_data/small.bam"));
  BOOST_CHECK(r.Close("test_data/small.cram"));
  
  // problematic here FIXME
  //SeqLib::BamReader r2;
  //BOOST_CHECK(r2.Open("tmp_out_poly.bam"));
  // should return false, no index
  //BOOST_CHECK(!r2.SetRegion(SeqLib::GenomicRegion(r.Header().Name2ID("X"),1001000, 1001100)));

}


BOOST_AUTO_TEST_CASE( plot_test ) {

  SeqLib::BamReader r;
  r.Open("test_data/small.bam");
  
  // should return false on empty region
  BOOST_CHECK(!r.SetMultipleRegions(SeqLib::GRC()));

  SeqLib::GenomicRegion gr("X:1,002,942-1,003,294", r.Header());
  r.SetRegion(gr);

  SeqLib::SeqPlot s;

  s.SetView(gr);

  SeqLib::BamRecord rec;
  SeqLib::BamRecordVector brv;
  while(r.GetNextRecord(rec))
    if (!rec.CountNBases() && rec.MappedFlag())
      brv.push_back(rec);

  s.SetPadding(20);

  std::cout << s.PlotAlignmentRecords(brv);

}


// CURRENTLY DOES NOT WORK
// need to find how to do reset
// BOOST_AUTO_TEST_CASE ( reset_works ) {

//   SeqLib::BamReader r;
//   r.Open("test_data/small.bam");
//   //r.Open("test_data/small.cram");

//   SeqLib::BamRecord rec1, rec2;
//   r.GetNextRecord(rec1);
//   r.Reset();
//   std::cerr << " AFTER RESET " << std::endl;
//   std::cerr << r.GetNextRecord(rec2) << std::endl;

//   BOOST_CHECK_EQUAL(rec1.Qname(), rec2.Qname());
//   }

BOOST_AUTO_TEST_CASE (json_parse) {

  SeqLib::BamReader r;
  r.Open("test_data/small.bam");
  ReadFilterCollection rfc(JSON1, r.Header());

  ReadFilter rf;
  SeqLib::GRC g(VCFFILE, r.Header());
  rf.addRegions(g);
  AbstractRule ar;
  ar.isize = Range(10,100, false);
  rf.SetMateLinked(true);
  rf.AddRule(ar);
  rfc.AddReadFilter(rf);

  std::cout << rfc << std::endl;

  SeqLib::BamRecord rec;
  size_t count = 0;
  int start, end;
  while(r.GetNextRecord(rec) && ++count < 10) {
    rec.QualityTrimmedSequence(4, start, end); // phred trim first
    rfc.isValid(rec);
  }

  // empty
  ReadFilterCollection rfc2("", r.Header());
  
}


BOOST_AUTO_TEST_CASE ( ref_genome ) {

  //SeqLib::RefGenome r("test_data/test_ref.fa");
  SeqLib::RefGenome r;
  r.LoadIndex("test_data/test_ref.fa");

  BOOST_CHECK(!r.IsEmpty());

  std::string out = r.QueryRegion("ref1", 0, 5);
  BOOST_CHECK_EQUAL(out, "ATCGAC");

  BOOST_CHECK_THROW(r.QueryRegion("ref1", 5,4), std::invalid_argument);
  BOOST_CHECK_THROW(r.QueryRegion("ref1", -1,4), std::invalid_argument);

  SeqLib::RefGenome r2;
  BOOST_CHECK_THROW(r2.QueryRegion("ref1",1,2), std::invalid_argument);

  // reload
  r2.LoadIndex("test_data/test_ref.fa");
}

BOOST_AUTO_TEST_CASE ( set_cigar ) {

  SeqLib::BamReader rr;
  rr.Open(SBAM); 
  SeqLib::BamRecord rec;
  size_t count = 0;
  while (rr.GetNextRecord(rec) && ++count < 10) {
      SeqLib::Cigar c;
      c.add(SeqLib::CigarField('M', 70));
      c.add(SeqLib::CigarField('I', 80));
      c.add(SeqLib::CigarField('M',1));
      rec.SetCigar(c);
      std::cerr << rec << std::endl;
  }


}
