#include "SeqLib/ReadFilter.h"

#include <cassert>
#include "htslib/htslib/khash.h"

//#define QNAME "D0EN0ACXX111207:7:2306:6903:136511"
//#define QFLAG -1

//#define DEBUG_MINI 1

#ifdef QNAME
#define DEBUGIV(msg, read)				\
  if (read.Qname() == QNAME && (read.AlignmentFlag() == QFLAG || QFLAG == -1)) { std::cerr << (msg) << " read " << r << std::endl; }
#else
#define DEBUGIV(msg, read)
#endif

namespace SeqLib {

  namespace Filter {
  // return if this rule accepts all reads
  bool AbstractRule::isEvery() const {
    return read_group.empty() && ins.isEvery() && del.isEvery() && isize.isEvery() && 
      mapq.isEvery() && len.isEvery() && clip.isEvery() && nm.isEvery() && 
      nbases.isEvery() && fr.isEvery() && 
      (subsam_frac >= 1) && xp.isEvery() 
#ifdef HAVE_C11
      && !aho.count
#endif
      ;
  }

// define what is a valid condition
/*    static const StringSet valid = 
  { 
  "duplicate", "supplementary", "qcfail", "hardclip", "fwd_strand",
  "rev_strand", "mate_fwd_strand", "mate_rev_strand", "mapped",
  "mate_mapped", "isize","clip", "length","nm",
  "mapq", "all", "ff", "xp","fr","rr","rf",
  "ic", "discordant","motif","nbases","!motif","allflag", "!allflag", "anyflag", "!anyflag",
  "ins","del",  "subsample", "rg"
};

    static const StringSet allowed_region_annots = 
    { "region","pad", "matelink", "exclude", "rules"};

    static const StringSet allowed_flag_annots = 
    {"duplicate", "supplementary", "qcfail", "hardclip", 
     "fwd_strand", "rev_strand", "mate_fwd_strand", "mate_rev_strand",
     "mapped", "mate_mapped", "ff", "fr", "rr", "rf", "ic"};
*/
bool ReadFilter::isValid(const BamRecord &r) {

  // empty default is pass
  if (!m_abstract_rules.size())
    return true;

  for (std::vector<AbstractRule>::iterator it = m_abstract_rules.begin(); 
       it != m_abstract_rules.end(); ++it) 
    if (it->isValid(r)) {
      ++it->m_count; //update this rule counter
      ++m_count;
       return true; // it is includable in at least one. 
    }
      
  return false;

}

  int FlagRule::parse_json_int(const Json::Value& v) {
    
    if (v.asInt())
      return v.asInt();
    else
      throw std::invalid_argument("Failed to parse int tag in JSON");
    return 0;
  }

  bool convert_to_bool(const Json::Value& value, const std::string& name) {

    Json::Value null(Json::nullValue);
    Json::Value v = value.get(name, null);
    if (v != null) {
      if (v.asBool())
	return true;
      else if (!v.asBool())
	return false;
    }

    return false;
    
  }

// check whether a BamAlignment (or optionally it's mate) is overlapping the regions
// contained in these rules
bool ReadFilter::isReadOverlappingRegion(const BamRecord &r) const {

  // if this is a whole genome rule, it overlaps
  if (!m_grv.size()) 
    return true;

  if (m_grv.CountOverlaps(GenomicRegion(r.ChrID(), r.Position(), r.PositionEnd())))
    return true;
  
  if (!m_applies_to_mate)
    return false;
  if (m_grv.CountOverlaps(GenomicRegion(r.MateChrID(), r.MatePosition(), r.MatePosition() + r.Length())))
    return true;

  return false;
}

// checks which rule a read applies to (using the hiearchy stored in m_regions).
// if a read does not satisfy a rule it is excluded.
  bool ReadFilterCollection::isValid(const BamRecord &r) {

  ++m_count_seen;

  DEBUGIV(r, "starting RFC isValid")

  if (m_regions.size() == 0)
    return true;

  DEBUGIV(r, "starting RFC isValid with non-empty regions")
  
    bool is_valid = false;
    bool exclude_hit = false;

    for (std::vector<ReadFilter>::iterator it = m_regions.begin(); it != m_regions.end(); ++it) {
      
      // only check read validity if it overlaps region
      if (!it->isReadOverlappingRegion(r)) 
	continue;
    
    // check the region with all its rules
      if (it->isValid(r)) {
     
      // if this is excluder region, exclude read
      if (it->excluder)
	exclude_hit = true;
      
      // in case we do fall through, track that we passed here
      is_valid = true;

    }
  }

    // found a hit in a rule
    if (is_valid && !exclude_hit) {
      ++m_count;
      return true;
    }
    
    return false;
}

  void ReadFilter::AddRule(const AbstractRule& ar) {
    m_abstract_rules.push_back(ar);
  }

  // constructor to make a ReadFilterCollection from a rules file.
  // This will reduce each individual BED file and make the 
  // GenomicIntervalTreeMap
  ReadFilterCollection::ReadFilterCollection(const std::string& script, const BamHeader& hdr) : m_count(0), m_count_seen(0) {

    // if is a file, read into a string
    std::ifstream iscript(script.c_str());
    std::stringstream ss;
    std::string tscript = script;
    if (iscript.is_open()) {
      char c = iscript.get();
      while (iscript.good()) {
	ss << c;
	c = iscript.get();
      }
      tscript = ss.str();;
    }

    // set up JsonCPP reader and attempt to parse script
    Json::Value root;
    Json::Reader reader;
    if ( !reader.parse(tscript, root)) {
      if (script.empty()) {
	std::cerr << "JSON script is empty. Setting default to filter all reads" << std::endl;
	return;
      }
      throw std::invalid_argument("ERROR: failed to parse JSON script");
    }

    Json::Value null(Json::nullValue);
    
    int level = 1;

    // assign the global rule if there is one
    // remove from the rest of the rules
    Json::Value glob = root.removeMember("global");
    if (!glob.isNull()) 
      rule_all.parseJson(glob);
    
    // iterator over regions
    for (Json::ValueConstIterator regions = root.begin(); regions != root.end(); ++regions) {

      ReadFilter mr;
      
      // check if mate applies
      mr.m_applies_to_mate = convert_to_bool(*regions, "matelink");

      // check for region padding
      int pad = regions->get("pad", 0).asInt();
      
      // set the region
      std::string reg;
      Json::Value v  = regions->get("region", null);
      if (v != null) {
	reg = v.asString();
	mr.id = mr.id + reg;
      }
      
      // actually parse the region
      if (reg == "WG" || reg.empty())
	mr.m_grv.clear(); // ensure it is whole-genome
      else {
	GRC regr(reg, hdr);
	regr.Pad(pad);
	mr.setRegions(regr);
      }
	// debug mr.setRegionFromFile(reg, hdr);
      
      // check if its excluder region
      mr.excluder = false; // default is no exclude
      v = regions->get("exclude", null);
      if (v != null) {
	mr.excluder = v.asBool();
	if (mr.excluder)
	  mr.id = mr.id + "_exclude";
      }
      
      // set the rules
      v = regions->get("rules", null);
      if (!v.size()) {
	//std::cerr << " !!!! RULES size is zero. !!!! " << std::endl;
	//exit(EXIT_FAILURE);
      }

      // loop through the rules
      for (Json::ValueIterator vv = v.begin(); vv != v.end(); ++vv) {
	if (*vv != null) {
	  AbstractRule ar = rule_all; // always start with the global rule
	  ar.parseJson(*vv);
	  // add the rule to the region
	  mr.m_abstract_rules.push_back(ar);
	}
      }
      
      // check that the regions have at least one rule
      // if it it doesn't, give it the global WG all
      if (!mr.m_abstract_rules.size())
	mr.m_abstract_rules.push_back(rule_all);
      
      mr.id = tostring(level);

      m_regions.push_back(mr);
      
    }
    
    // check that there is at least one non-excluder region. 
    // if not, give global includer
    bool has_includer = false;
    for (std::vector<ReadFilter>::const_iterator kk = m_regions.begin(); kk != m_regions.end(); ++kk) 
      if (!kk->excluder)
	has_includer = true;
    if (!has_includer) {
      ReadFilter mr;
      mr.m_abstract_rules.push_back(rule_all);
      mr.id = "WG_includer";
      m_regions.push_back(mr);
    }
    
  }
  
  void ReadFilter::setRegions(const GRC& g) {
    m_grv = g;
    m_grv.CreateTreeMap();
  }

  void ReadFilter::addRegions(const GRC& g) {
    m_grv.Concat(g);
    m_grv.MergeOverlappingIntervals();
    m_grv.CreateTreeMap();
  }


  // print the ReadFilterCollection
  std::ostream& operator<<(std::ostream &out, const ReadFilterCollection &mr) {
    
    out << "----------ReadFilterCollection-------------" << std::endl;

    for (std::vector<ReadFilter>::const_iterator it = mr.m_regions.begin(); it != mr.m_regions.end(); ++it) 
      out << *it << std::endl;
    out << "------------------------------------------";
    return out;
    
  }

// print a ReadFilter information
std::ostream& operator<<(std::ostream &out, const ReadFilter &mr) {
  
  std::string file_print = !mr.m_grv.size() ? "WHOLE GENOME" : mr.m_region_file;
  out << (mr.excluder ? "--Exclude Region: " : "--Include Region: ") << file_print;
  if (mr.m_grv.size()) {
    //out << " --Size: " << AddCommas<int>(mr.m_width); 
    out << " Matelink: " << (mr.m_applies_to_mate ? "ON" : "OFF");
    if (mr.m_grv.size() == 1)
      out << " Region : " << mr.m_grv[0] << std::endl;
    else
      out << " " << mr.m_grv.size() << " regions" << std::endl;      
  } else {
    out << std::endl;
  }

  for (std::vector<AbstractRule>::const_iterator it = mr.m_abstract_rules.begin(); it != mr.m_abstract_rules.end(); ++it) 
    out << *it << std::endl;
  return out;
}

  void ReadFilterCollection::AddReadFilter(const ReadFilter& rf) {
    m_regions.push_back(rf);
  }

  ReadFilter::~ReadFilter() {}

  bool Flag::parseJson(const Json::Value& value, const std::string& name) {

    if (value.isMember(name.c_str())) {
      convert_to_bool(value, name) ? setOn() : setOff();
      return true;
    }
    
    return false;

  }

  void FlagRule::parseJson(const Json::Value& value) {

    Json::Value null(Json::nullValue);
    if (value.isMember("allflag"))
      setAllOnFlag(parse_json_int(value.get("allflag", null)));
    if (value.isMember("!allflag"))
      setAllOffFlag(parse_json_int(value.get("!allflag", null)));
    if (value.isMember("anyflag"))
      setAnyOnFlag(parse_json_int(value.get("anyflag", null)));
    if (value.isMember("!anyflag"))
      setAnyOffFlag(m_any_off_flag = parse_json_int(value.get("!anyflag", null)));
    
    // have to set the every if find flag so that rule knows it cant skip checking
    if (dup.parseJson(value, "duplicate")) every = false;
    if (supp.parseJson(value, "supplementary")) every = false;
    if (qcfail.parseJson(value, "qcfail")) every = false;
    if (hardclip.parseJson(value, "hardclip")) every = false;
    if (fwd_strand.parseJson(value, "fwd_strand")) every = false;
    if (mate_rev_strand.parseJson(value, "mate_rev")) every = false;
    if (mate_fwd_strand.parseJson(value, "mate_fwd")) every = false;
    if (mate_mapped.parseJson(value, "mate_mapped")) every = false;
    if (mapped.parseJson(value, "mapped")) every = false;
    if (ff.parseJson(value, "ff")) every = false;
    if (fr.parseJson(value, "fr")) every = false;
    if (rf.parseJson(value, "rf")) every = false;
    if (rr.parseJson(value, "rr")) every = false;
    if (ic.parseJson(value, "ic")) every = false;

  }
  
  void Range::parseJson(const Json::Value& value, const std::string& name) {
    Json::Value null(Json::nullValue);
    Json::Value v = value.get(name, null);

    if (v != null) {
      if (v.size() > 2) {
	std::cerr << " ERROR. Not expecting array size " << v.size() << " for Range " << name << std::endl;
      } else {
	m_every = false;
	m_inverted = false;

	if (v.isArray()) {
	  m_min = v[0].asInt();
	  m_max = v[1].asInt();
	} else if (v.isInt()) {
	  m_min = v.asInt();
	  m_max = INT_MAX;
	} else if (v.isBool()) {
	  m_min = v.asBool() ? 1 : INT_MAX; // if true, [1,MAX], if false [MAX,1] (not 1-MAX)
	  m_max = v.asBool() ? INT_MAX : 1;
	} else {
	  throw std::invalid_argument("Unexpected type for range flag: " + name);
	}

	if (m_min > m_max) {
	  m_inverted = true;
	  std::swap(m_min, m_max); // make min always lower
	}
      }
	
    }
  }

  void AbstractRule::parseJson(const Json::Value& value) {

    // parse read group
    const std::string rg = "rg";
    if (value.isMember(rg.c_str())) {
      Json::Value null(Json::nullValue);
      Json::Value v = value.get(rg, null);
      assert(v != null);
      read_group = v.asString();
    }
      
    // set the ID
    std::vector<std::string> mn = value.getMemberNames();
    for (std::vector<std::string>::const_iterator i = mn.begin(); i != mn.end(); ++i)
      id += *i + ";";

    // not necessary, not c++98 compatible
    //if (id.length())
    //  id.pop_back();

    // parse the flags
    fr.parseJson(value);
    
    isize.parseJson(value, "isize");
    mapq.parseJson(value, "mapq");
    len.parseJson(value, "length");
    clip.parseJson(value, "clip");
    nbases.parseJson(value, "nbases");
    ins.parseJson(value, "ins");
    del.parseJson(value, "del");
    nm.parseJson(value, "nm");
    xp.parseJson(value, "xp");
    
    // parse the subsample data
    parseSubLine(value);

    // parse the motif line
    parseSeqLine(value);
    
  }


    // main function for determining if a read is valid
    bool AbstractRule::isValid(const BamRecord &r) {
    
      DEBUGIV(r, "starting AR:isValid")

    // check if its keep all or none
    if (isEvery())
      return true;
    
    // check if it is a subsample
    if (subsam_frac < 1) {
      uint32_t k = __ac_Wang_hash(__ac_X31_hash_string(r.QnameChar()) ^ subsam_seed);
      if ((double)(k&0xffffff) / 0x1000000 >= subsam_frac) 
	return false;
    }
    
    // check if is discordant
    bool isize_pass = isize.isValid(r.FullInsertSize());

    DEBUGIV("isize_pass: " + isize_pass, r)
    
    if (!isize_pass) {
      return false;
    }
    
    // check for valid read name 
    if (!read_group.empty()) {
      std::string RG = r.ParseReadGroup();
      if (!RG.empty() && RG != read_group)
	return false;
    }

    // check for valid mapping quality
    if (!mapq.isEvery())
      if (!mapq.isValid(r.MapQuality())) 
	return false;

    DEBUGIV(r, "mapq pass")
    
    // check for valid flags
    if (!fr.isValid(r)) 
      return false;

    DEBUGIV(r, "flag pass")
    
    // check the CIGAR
    if (!ins.isEvery() || !del.isEvery()) {
      if (!ins.isValid(r.MaxInsertionBases()))
	return false;
      if (!del.isValid(r.MaxDeletionBases()))
	return false;
    }

    DEBUGIV(r, "cigar pass")
      
    // get the sequence as trimmed
    std::string tseq = r.QualitySequence(); //AddZTag("GV", r.Sequence().substr(startpoint, new_len));
    
#ifdef HAVE_C11
    // check for aho corasick motif match
    if (aho.count) {
      if (!aho.QueryText(tseq))
      return false;
      DEBUGIV(r, "aho pass")
    }
#endif    

    // check for valid NM
    if (!nm.isEvery()) {
      int32_t nm_val = 0;
      r.GetIntTag("NM", nm_val);
      if (!nm.isValid(nm_val))
	return false;
      DEBUGIV(r, "NM pass")
    }
    
    // check the N bases
    if (!nbases.isEvery()) {
      size_t n = r.CountNBases();
      if (!nbases.isValid(n))
	return false;
      DEBUGIV(r, "N bases pass")
    }

    // check for valid length
    if (!len.isValid(tseq.length())) {
      return false;
      DEBUGIV(r, "len pass")
    }

    // check for valid clip
    int new_clipnum = r.NumClip() - (r.Length() - tseq.length()); // get clips, minus amount trimmed off
    if (!clip.isValid(new_clipnum)) {
      return false;
      DEBUGIV(r, "clip pass with clip size " + tostring(new_clipnum))
    }

    // check for secondary alignments
    if (!xp.isEvery()) {
      if (!xp.isValid(r.CountBWASecondaryAlignments())) {
	return false;
      }
      DEBUGIV(r, "XP pass")
    }

    DEBUGIV(r, "**** READ ACCEPTED IN AR:ISVALID")
    return true;
  }
  
  bool FlagRule::isValid(const BamRecord &r) {
    
    DEBUGIV(r, "flagrule start")

    if (isEvery())
      return true;

    // if have on or off flag, use that
    // 0001100 - all flag
    // 0101000 - flag
    // -------
    // 0001000 - should fail all flag. Should pass any flag
    if (m_all_on_flag && !( (r.AlignmentFlag() & m_all_on_flag)  == m_all_on_flag) ) // if all on, pass
      return false;
    if (m_all_off_flag && ( (r.AlignmentFlag() & m_all_off_flag) == m_all_off_flag) ) // if all on, fail
      return false;

    // if have on or off flag, use that
    if (m_any_on_flag && !(r.AlignmentFlag() & m_any_on_flag) ) // if ANY on, pass
      return false;
    if (m_any_off_flag && (r.AlignmentFlag() & m_any_off_flag)) // if ANY on, fail
      return false;

    DEBUGIV(r, "FlagRule::isValid checking named flags")

    if (!dup.isNA()) 
      if ((dup.isOff() && r.DuplicateFlag()) || (dup.isOn() && !r.DuplicateFlag()))
      return false;
  if (!supp.isNA()) 
    if ((supp.isOff() && r.SecondaryFlag()) || (supp.isOn() && !r.SecondaryFlag()))
      return false;
  if (!qcfail.isNA())
    if ((qcfail.isOff() && r.QCFailFlag()) || (qcfail.isOn() && !r.QCFailFlag()))
      return false;
  if (!mapped.isNA())
    if ( (mapped.isOff() && r.MappedFlag()) || (mapped.isOn() && !r.MappedFlag()))
      return false;
  if (!mate_mapped.isNA())
    if ( (mate_mapped.isOff() && r.MateMappedFlag()) || (mate_mapped.isOn() && !r.MateMappedFlag()) )
      return false;

  // check for hard clips
  if (!hardclip.isNA())  {// check that we want to chuck hard clip
    if (r.CigarSize() > 1) {
      bool ishclipped = r.NumHardClip() > 0;
      if ( (ishclipped && hardclip.isOff()) || (!ishclipped && hardclip.isOn()) )
	return false;
    }
  }

  // check for orientation
  // check first if we need to even look for orientation
  bool ocheck = !ff.isNA() || !fr.isNA() || !rf.isNA() || !rr.isNA() || !ic.isNA();

  // now its an orientation pair check. If not both mapped, chuck
  if (!r.PairMappedFlag() && ocheck)
    return false;

  if ( ocheck ) {

    //    bool first = r.Position() < r.MatePosition();
    //bool bfr = (first && (!r.ReverseFlag() && r.MateReverseFlag())) || (!first &&  r.ReverseFlag() && !r.MateReverseFlag());
    //bool brr = r.ReverseFlag() && r.MateReverseFlag();
    //bool brf = (first &&  (r.ReverseFlag() && !r.MateReverseFlag())) || (!first && !r.ReverseFlag() &&  r.MateReverseFlag());
    //bool bff = !r.ReverseFlag() && !r.MateReverseFlag();
      
    bool bic = r.Interchromosomal();

    int PO = r.PairOrientation();

    // its FR and it CANT be FR (off) or its !FR and it MUST be FR (ON)
    // orienation not defined for inter-chrom, so exclude these with !ic
    if (!bic) { // PROCEED IF INTRA-CHROMOSOMAL
      //if ( (bfr && fr.isOff()) || (!bfr && fr.isOn())) 
      if ( (PO == FRORIENTATION && fr.isOff()) || (PO != FRORIENTATION && fr.isOn())) 
	return false;
      //if ( (brr && rr.isOff()) || (!brr && rr.isOn())) 
      if ( (PO == RRORIENTATION && rr.isOff()) || (PO != RRORIENTATION && rr.isOn())) 
	return false;
      //if ( (brf && rf.isOff()) || (!brf && rf.isOn())) 
      if ( (PO == RFORIENTATION && rf.isOff()) || (PO != RFORIENTATION&& rf.isOn())) 
	return false;
      //if ( (bff && ff.isOff()) || (!bff && ff.isOn())) 
      if ( (PO == FFORIENTATION && ff.isOff()) || (PO != FFORIENTATION && ff.isOn())) 
	return false;
    }
    if ( (bic && ic.isOff()) || (!bic && ic.isOn()))
      return false;
      
  }

  return true;
  
}

// define how to print
std::ostream& operator<<(std::ostream &out, const AbstractRule &ar) {

  out << "  Rule: ";
  if (ar.isEvery()) {
    out << "  ALL";
  } else {
    if (!ar.read_group.empty())
      out << "Read Group: " << ar.read_group << " -- ";
    if (!ar.isize.isEvery())
      out << "isize:" << ar.isize << " -- " ;
    if (!ar.mapq.isEvery())
      out << "mapq:" << ar.mapq << " -- " ;
    if (!ar.len.isEvery())
      out << "length:" << ar.len << " -- ";
    if (!ar.clip.isEvery())
      out << "clip:" << ar.clip << " -- ";
    if (!ar.nm.isEvery())
      out << "nm:" << ar.nm << " -- ";
    if (!ar.xp.isEvery())
      out << "xp:" << ar.xp << " -- ";
    if (!ar.nbases.isEvery())
      out << "nbases:" << ar.nbases << " -- ";
    if (!ar.ins.isEvery())
      out << "ins:" << ar.ins << " -- ";
    if (!ar.del.isEvery())
      out << "del:" << ar.del << " -- ";
    if (ar.subsam_frac < 1)
      out << "sub:" << ar.subsam_frac << " -- ";
#ifdef HAVE_C11
    if (ar.aho.count)
      out << "motif: " << ar.aho.file << " -- ";
#endif
    out << ar.fr;
  }
  return out;
}

// define how to print
std::ostream& operator<<(std::ostream &out, const FlagRule &fr) {

  if (fr.isEvery()) {
    out << "Flag: ALL";
    return out;
  } 

  std::string keep = "Flag ON: ";
  std::string remo = "Flag OFF: ";

  if (fr.m_all_on_flag)
    keep += "[(all)" + tostring(fr.m_all_on_flag) + "],";
  if (fr.m_all_off_flag)
    remo += "[(all)" + tostring(fr.m_all_off_flag) + "],";

  if (fr.m_any_on_flag)
    keep += "[(any)" + tostring(fr.m_any_on_flag) + "],";
  if (fr.m_any_off_flag)
    remo += "[(any)" + tostring(fr.m_any_off_flag) + "],";

  if (fr.dup.isOff())
    remo += "duplicate,";
  if (fr.dup.isOn())
    keep += "duplicate,";

  if (fr.supp.isOff())
    remo += "supplementary,";
  if (fr.supp.isOn())
    keep += "supplementary,";

  if (fr.qcfail.isOff())
    remo += "qcfail,";
  if (fr.qcfail.isOn())
    keep += "qcfail,";

  if (fr.hardclip.isOff())
    remo += "hardclip,";
  if (fr.hardclip.isOn())
    keep += "hardclip,";

  if (fr.paired.isOff())
    remo += "paired,";
  if (fr.paired.isOn())
    keep += "paired,";



  if (fr.ic.isOff())
    remo += "ic,";
  if (fr.ic.isOn())
    keep += "ic,";

  if (fr.ff.isOff())
    remo += "ff,";
  if (fr.ff.isOn())
    keep += "ff,";

  if (fr.fr.isOff())
    remo += "fr,";
  if (fr.fr.isOn())
    keep += "fr,";

  if (fr.rr.isOff())
    remo += "rr,";
  if (fr.rr.isOn())
    keep += "rr,";

  if (fr.rf.isOff())
    remo += "rf,";
  if (fr.rf.isOn())
    keep += "rf,";

  if (fr.mapped.isOff())
    remo += "mapped,";
  if (fr.mapped.isOn())
    keep += "mapped,";

  if (fr.mate_mapped.isOff())
    remo += "mate_mapped,";
  if (fr.mate_mapped.isOn())
    keep += "mate_mapped,";

  keep = keep.length() > 10 ? keep.substr(0, keep.length() - 1) : ""; // remove trailing comment
  remo = remo.length() > 10 ? remo.substr(0, remo.length() - 1) : ""; // remove trailing comment
  
  if (!keep.empty() && !remo.empty())
    out << keep << " -- " << remo;
  else if (!keep.empty())
    out << keep;
  else
    out << remo;

  return out;
}

// define how to print
std::ostream& operator<<(std::ostream &out, const Range &r) {
  if (r.isEvery())
    out << "ALL";
  else
    out << (r.m_inverted ? "NOT " : "") << "[" << r.m_min << "," << (r.m_max == INT_MAX ? "MAX" : tostring(r.m_max))  << "]";
  return out;
}

  void AbstractRule::parseSeqLine(const Json::Value& value) {
    
#ifdef HAVE_C11
    bool i = false; // invert motif?
#endif
    std::string motif_file;
    Json::Value null(Json::nullValue);
    if (value.get("motif", null) != null) 
      motif_file = value.get("motif", null).asString();
    else if (value.get("!motif", null) != null) {
      motif_file = value.get("!motif", null).asString();
#ifdef HAVE_C11
      i = true;
#endif
    }
    else
      return;
#ifdef HAVE_C11
    addMotifRule(motif_file, i);
#else
    if (!motif_file.empty())
      std::cerr << "WARNING: Need to compile with C++11 for Aho-Corasick matching" << std::endl;
#endif

  return;

  }

#ifdef HAVE_C11
  void AbstractRule::addMotifRule(const std::string& f, bool inverted) {
    std::cerr << "...making the AhoCorasick trie from " << f << std::endl;
    aho.TrieFromFile(f);
    std::cerr << "...finished making AhoCorasick trie with " << AddCommas(aho.count) << " motifs" << std::endl;
    aho.inv = inverted;
  }

  void AhoCorasick::TrieFromFile(const std::string& f) {

    file = f;

    // open the sequence file
    std::ifstream iss(f.c_str());
    if (!iss || !read_access_test(f)) 
      throw std::runtime_error("AhoCorasick::TrieFromFile - Cannot read file: " + f);
    
    // make the Aho-Corasick trie
    std::string pat;
    while (getline(iss, pat, '\n')) {
      ++count;
      AddMotif(pat);
    }
  }
#endif
  
  void AbstractRule::parseSubLine(const Json::Value& value) {
    Json::Value null(Json::nullValue);
    if (value.get("subsample", null) != null) 
      subsam_frac = value.get("subample", null).asDouble();
  }
  
GRC ReadFilterCollection::getAllRegions() const
{
  GRC out;

  for (std::vector<ReadFilter>::const_iterator i = m_regions.begin(); i != m_regions.end(); ++i)
    out.Concat(i->m_grv);

  return out;
}
    
#ifdef HAVE_C11
    int AhoCorasick::QueryText(const std::string& t) const {
      auto matches = aho_trie->parse_text(t);
      return matches.size();
      return 0;
    }
#endif

  }
}
