#ifndef UTIL_HH
#define UTIL_HH 1

/****************************************************************************
** Util.hh
**
** Routines for IO and DNA sequence analysis
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

#include <string>
#include <set>
#include <sstream>


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
std::string itos(int i);
std::string dtos(double d);
bool isDNA(char b);
char rrc(char b);
std::string rc_str(const std::string & str);
bool Fasta_Read(FILE * fp, std::string & s, std::string & hdr);
bool isNseq(const std::string & seq);
bool isRepeat(const std::string & seq, int K);
bool isAlmostRepeat(const std::string & seq, int K, int max);
bool kMismatch(size_t s, size_t e, const std::string & t, size_t start, int max);
//void findTandems(const std::string & seq, const string & tag);

#endif
