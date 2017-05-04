#include "SeqLib/SeqPlot.h"

namespace SeqLib {

std::string SeqPlot::PlotAlignmentRecords(const BamRecordVector& brv) const {

  PlottedReadVector plot_vec;

  for (BamRecordVector::const_iterator i = brv.begin(); i != brv.end(); ++i) {
    
    // get the position in the view window
    if (i->ChrID() != m_view.chr)
      continue;

    int pos = i->Position() - m_view.pos1;
    if (pos < 0)
      continue;

    if (i->PositionEnd() > m_view.pos2)
      continue;

    // plot with gaps
    std::string tseq = i->Sequence();
    std::string gapped_seq;

    size_t p = i->AlignmentPosition(); // move along on sequence, starting at first non-clipped base
    Cigar cc = i->GetCigar();
    for (Cigar::const_iterator c = cc.begin(); c != cc.end(); ++c) {
      if (c->Type() == 'M') { // 
	assert(p + c->Length() <= tseq.length());
	gapped_seq += tseq.substr(p, c->Length());
      } else if (c->Type() == 'D') {
	gapped_seq += std::string(c->Length(), '-');
      }

      if (c->Type() == 'I' || c->Type() == 'M')
	p += c->Length();
    }

    std::stringstream msg;
    msg << i->Qname() << ">>>" << (i->ChrID() + 1) << ":" << i->Position();
      
    // add to the read plot
    plot_vec.push_back(PlottedRead(pos, gapped_seq, msg.str()));
    
  }

  // sort them
  std::sort(plot_vec.begin(), plot_vec.end());

  // make a list of lines
  PlottedReadLineVector line_vec;

  // plot the reads from the ReadPlot vector
  for (PlottedReadVector::iterator i = plot_vec.begin(); i != plot_vec.end(); ++i) {
    bool found = false;
    for (PlottedReadLineVector::iterator j = line_vec.begin(); j != line_vec.end(); ++j) {
      if (j->readFits(*i)) { // it fits here
	j->addRead(&(*i));
	found = true;
	break;
      }
    }
    if (!found) { // didn't fit anywhere, so make a new line
      PlottedReadLine prl;
      prl.pad = m_pad;
      prl.contig_len = m_view.Width(); //ac.getSequence().length();
      prl.addRead(&(*i));
      line_vec.push_back(prl);
    }
  }
  
  std::stringstream ss;
  
  // plot the lines. Add contig identifier to each
  for (PlottedReadLineVector::const_iterator i = line_vec.begin(); i != line_vec.end(); ++i) 
    ss << (*i) << std::endl;

  return ss.str();

  
}


}
