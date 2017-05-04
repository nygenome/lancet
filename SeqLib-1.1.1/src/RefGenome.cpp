#include "SeqLib/RefGenome.h"

#include <stdexcept>
#include "SeqLib/SeqLibUtils.h"

namespace SeqLib {

  bool RefGenome::LoadIndex(const std::string& file) {

    // clear the old one
    if (index)  
      fai_destroy(index);
    
    index = NULL;
    
    // check that its readable
    if (!read_access_test(file)) {
      return false;
      //throw std::invalid_argument("RefGenome: file not found - " + file);
    }
    
    // load it in
    index = fai_load(file.c_str());

    if (!index)
      return false;

    return true;

  }
  
  std::string RefGenome::QueryRegion(const std::string& chr_name, int32_t p1, int32_t p2) const {
    
    // check that we have a loaded index
    if (!index) 
      throw std::invalid_argument("RefGenome::queryRegion index not loaded");

    // check input is OK
    if (p1 > p2)
      throw std::invalid_argument("RefGenome::queryRegion p1 must be <= p2");
    if (p1 < 0)
      throw std::invalid_argument("RefGenome::queryRegion p1 must be >= 0");

    int len;
    char * f = faidx_fetch_seq(index, const_cast<char*>(chr_name.c_str()), p1, p2, &len);

    if (!f)
      throw std::invalid_argument("RefGenome::queryRegion - Could not find valid sequence");

    std::string out(f);

    free(f);

    if (out.empty())
      throw std::invalid_argument("RefGenome::queryRegion - Returning empty query on " + chr_name + ":" + tostring(p1) + "-" + tostring(p2));

    return (out);

  }

}
