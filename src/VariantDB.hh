#ifndef VARIANTDB_HH
#define VARIANTDB_HH 1

/******************************************************************
** VariantDB.hh
**
** Class for storing multiple variants (DB-style)
**
**  Authors: Giuseppe Narzisi
**    Date: December 21, 2015
**
*******************************************************************/

#include <map>
#include <string>
#include <iostream>
#include "util.hh"
#include "Variant.hh"
#include "sha256.hh"

using namespace std;

class VariantDB_t
{
public:

	map<string,Variant_t> DB;

	VariantDB_t() {}
	
	void addVar(Variant_t v);
	void printToVCF();
};

#endif
