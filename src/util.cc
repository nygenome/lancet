#include "util.hh"
#include <iostream>
#include <cstdio>

/******************************************************************
** Util.cc
**
** Routines for IO and DNA sequence analysis
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

using namespace std;


// xfopen
//////////////////////////////////////////////////////////////////////////

FILE * xfopen(const string & filename, const string & mode)
{
	if (filename == "-")
	{
		if (mode == "r") { return stdin; }
		else             { return stdout; }
	}

	FILE * fp = fopen(filename.c_str(), mode.c_str());

	if (!fp)
	{
		cerr << "Couldn't open " << filename << " for " << mode << endl;
		exit(1);
	}

	return fp;
}

// xfclose
//////////////////////////////////////////////////////////////////////////

void xfclose(FILE * fp)
{
	if (fp != stdin && fp != stdout)
	{
		fclose(fp);
	}
}


// isDNA
//////////////////////////////////////////////////////////////////////////

bool isDNA(char b)
{
	if (b == 'A' || b == 'a' ||
		b == 'C' || b == 'c' ||
		b == 'G' || b == 'g' ||
		b == 'T' || b == 't')
	{ 
		return true;
	}

	return false;
}

// rc
//////////////////////////////////////////////////////////////////////////

char rrc(char b)
{
	switch(b)
	{
		case 'A' : return 'T'; case 'a' : return 't';
		case 'C' : return 'G'; case 'c' : return 'g';
		case 'G' : return 'C'; case 'g' : return 'c';
		case 'T' : return 'A'; case 't' : return 'a';
		case 'N' : return 'N'; case 'n' : return 'n'; 
		
	}

	return 0;
}

// rc_str
//////////////////////////////////////////////////////////////////////////
string rc_str(const string & str)
{
	string retval;

	for (int i = str.length()-1; i >= 0; i--)
	{
		retval.push_back(rrc(str[i]));
	}

	return retval;
}

// isNseq
// returns true if the input sequence contains only Ns
//////////////////////////////////////////////////////////////////////////
bool isNseq(const string & seq)
{
	bool result = true;

	int end = seq.length();
	for (int i = 0; i < end; i++)
	{

		if ( (seq[i] != 'N') || (seq[i] != 'n') ) {
			result = false;
			break;
		}
	}
	return result;	
}

// isRepeat
// returns true if the input sequence contains any repeat of size K
// (multiple occurence of the same K-mer)
//////////////////////////////////////////////////////////////////////////
bool isRepeat(const string & seq, int K)
{
	bool result = false;
	
	set<string> mers;

	int end = seq.length() - K;
	for (int offset = 0; offset < end; offset++)
	{
		string s = seq.substr(offset,K);

	    if (mers.count(s)>0) { // s is in the set
			result = true;
			break;
		}
		else { // not in the set
			mers.insert(s);
		}
	}
	return result;	
}

bool isAlmostRepeat(const string & seq, int K, int max)
{
	bool result = false;
	
	int end = seq.length() - K;
	for (int offset = 0; offset < end; offset++)
	{
		//string s = seq.substr(offset,K);

		int start = offset;
		int end = offset + K;
		if(kMismatch(start,end,seq,offset+1,max)) {
			result = true;
			break;
		}
	}
	return result;	
}

bool kMismatch(size_t s, size_t e, const string & t, size_t start, int max) { // p==pattern, t==text
	bool flag;
	int count;
	size_t i=start;
	size_t L = e-s+1;
	//while(i<(t.size()-p.size()+1)) {
	while(i<(t.size()-L+1)) {
		flag = true;
		count = 0;
		size_t j=0;
		//while(j<p.size() && i+j<t.size()) {
		while(j<L && i+j<t.size()) {
			//if(t[i+j] != p[j]) {
			if(t[i+j] != t[s+j]) {
				count++;
				if(count>max) { flag = false; break; }
			}
			j++;
		}
		//if(flag && j==p.size()) { return true; }
		if(flag && j==L) { return true; }
		i++;
	}
	return false;
}


// Fasta_Read
//////////////////////////////////////////////////////////////////////////

bool Fasta_Read (FILE * fp, string & s, string & hdr)
{
	int  ch;

	s . erase ();
	hdr . erase ();

	// skip till next '>' if necessary
	while  ((ch = fgetc (fp)) != EOF && ch != '>')
		;

	if  (ch == EOF)
		return  false;

	// skip spaces if any
	while  ((ch = fgetc (fp)) != EOF && ch == ' ')
		;
	if  (ch == EOF)
		return  false;
	ungetc (ch, fp);

	// put rest of line into  hdr
	while  ((ch = fgetc (fp)) != EOF && ch != '\n')
		hdr . push_back (char (ch));

	// put everything up till next '>' into  s
	while  ((ch = fgetc (fp)) != EOF && ch != '>')
	{
		if  (! isspace (ch))
			s . push_back (char (ch));
	}

	if  (ch == '>')
		ungetc (ch, fp);

	return  true;
}

