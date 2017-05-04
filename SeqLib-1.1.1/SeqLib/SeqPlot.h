#ifndef SEQLIB_CONTIG_PLOT_H
#define SEQLIB_CONTIG_PLOT_H

#include "SeqLib/BamRecord.h"

namespace SeqLib {
  
  /** Object for creating ASCII alignment plots
   */
  class SeqPlot {

  public:
    
    /** Create an empty plotter object */
  SeqPlot() : m_pad(5) {}

    /** Plot aligned read by stacking them in an IGV-like view */
    std::string PlotAlignmentRecords(const BamRecordVector& brv) const;

    /** Set the view window 
     * @param g Window to view reads in. Reads that 
     * start or end outside of this window are not plotted
     */
    void SetView(const GenomicRegion& g) { m_view = g; }

    /** Set the padding between reads (default is 5) */
    void SetPadding(int p) { m_pad = p; }

  private: 

    // reads that align to the contig
    BamRecordVector m_reads;

    // view window
    GenomicRegion m_view;

    // padding when placing reads
    int m_pad;

  };

  /** A single plotted read */
struct PlottedRead {

  int pos;
  std::string seq;
  std::string info;
  
  PlottedRead(int p, const std::string& s, const std::string& i) : pos(p), seq(s), info(i) {}

  bool operator<(const PlottedRead& pr) const {
    return (pos < pr.pos);
  }

};

typedef std::vector<PlottedRead> PlottedReadVector;

/** A plotted line */
struct PlottedReadLine {

PlottedReadLine() : available(0), contig_len(0), pad(5) {}

  std::vector<PlottedRead*> read_vec;
  int available;
  int contig_len;
  
  int pad;

  void addRead(PlottedRead *r) {
    read_vec.push_back(r);
    available = r->pos + r->seq.length() + pad;
  }

  bool readFits(const PlottedRead &r) {
    return (r.pos >= available);
  }

  friend std::ostream& operator<<(std::ostream& out, const PlottedReadLine &r) {
    int last_loc = 0;
    for (std::vector<PlottedRead*>::const_iterator i = r.read_vec.begin(); i != r.read_vec.end(); ++i) {
      //    for (auto& i : r.read_vec) {
      assert((*i)->pos - last_loc >= 0);
      out << std::string((*i)->pos - last_loc, ' ') << (*i)->seq;
      last_loc = (*i)->pos + (*i)->seq.length();
    }
    int name_buff = r.contig_len - last_loc;
    assert(name_buff < 1e6);
    out << std::string(std::max(name_buff, 5), ' ');
    for (std::vector<PlottedRead*>::const_iterator i = r.read_vec.begin(); i != r.read_vec.end(); ++i) {
      //for (auto& i : r.read_vec) { // add the data
      out << (*i)->info << ",";
    }
    return out;
  }

};

typedef std::vector<PlottedReadLine> PlottedReadLineVector;


}



#endif
