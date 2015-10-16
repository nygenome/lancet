#ifndef UTIL_HH
#define UTIL_HH 1

/******************************************************************
** Util.hh
**
** Routines for IO and DNA sequence analysis
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <string>
#include <set>


//-- Include hash_map
#ifdef __GNUC__
#if __GNUC__ < 3
  #include <hash_map.h>
  namespace Sgi { using ::hash_map; }; // inherit globals
  #define HASHMAP std
#elif __GNUC__ == 3
  #include <ext/hash_map>
  #if __GNUC_MINOR__ == 0
    namespace Sgi = std;               // GCC 3.0
    #define HASHMAP std
  #else
    namespace Sgi = ::__gnu_cxx;       // GCC 3.1 and later
    #define HASHMAP __gnu_cxx
  #endif
#elif __GNUC__ > 3
  #include <ext/hash_map>
  namespace Sgi = ::__gnu_cxx;         // GCC 4.0 and later
  #define HASHMAP __gnu_cxx
#endif
#else      // ...  there are other compilers, right?
  namespace Sgi = std;
  #define HASHMAP std
#endif



FILE * xfopen(const std::string & filename, const std::string & mode);

void xfclose(FILE * fp);

bool isDNA(char b);

char rrc(char b);

std::string rc_str(const std::string & str);

bool Fasta_Read(FILE * fp, std::string & s, std::string & hdr);

bool isNseq(const std::string & seq);

bool isRepeat(const std::string & seq, int K);

bool isAlmostRepeat(const std::string & seq, int K, int max);

//bool kMismatch(const std::string & p, const std::string & t, int start, int max);
bool kMismatch(size_t s, size_t e, const std::string & t, size_t start, int max);

#endif
