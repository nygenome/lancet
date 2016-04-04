#include "util.hh"
#include <iostream>
#include <cstdio>

/****************************************************************************
** Util.cc
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

// integer to string conversion
//////////////////////////////////////////////////////////////////////////
string itos(int i) // convert int to string
{
    stringstream s;
    s << i;
    return s.str();
}

// double to string conversion
//////////////////////////////////////////////////////////////////////////
string dtos(double d) // convert int to string
{
    stringstream s;
    s << d;
    return s.str();
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

// reverse
// reverse the content of a string
//////////////////////////////////////////////////////////////
string reverse(const string & str) 
{
	string retval = str;
	
	int i=0;
	int j=retval.size()-1;
	while(i<j){
		char tmp = retval[i];
		retval[i] = retval[j];
		retval[j] = tmp; 
		i++;j--;
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

// returns true if all bases have phred quality >= Q 
//////////////////////////////////////////////////////////////////////////
bool seqAboveQual(string qv, int Q) 
{
	for ( string::const_iterator it=qv.begin(); it!=qv.end(); ++it) {
		if( *it < Q ) { return false; }
	}
	return true;
}

// parse MD string 
// extract SNVs locations and add them to map M
//////////////////////////////////////////////////////////////////////////
void parseMD(string & md, map<int,int> & M, int start) {
	// String for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*10
	// example: 6G4C20G1A5C5A1^C3A15G1G15
    map<int,int>::iterator mit;
	size_t p = md.find_first_of("ACGT^");
	size_t p_old = -1;
	size_t p2;
	string del;
	string num;
	int pos = start;
	while(p!=string::npos) {
		num = md.substr(p_old+1,p-(p_old+1));
		pos += atoi(num.c_str());
		//cerr << num << "|";
		if(md[p]=='^') { 
			p2 = md.find_first_not_of("ACGT^",p+1); // find next number
			del = md.substr(p,p2-p);
			pos += del.size();
			//cerr << del << "|";
			p = md.find_first_of("ACGT^",p2);
			p_old = p2-1;
		}
		else {
			pos++;
			mit = M.find(pos);
			if (mit != M.end()) { ((*mit).second)++; }
			else { M.insert(std::pair<int,int>(pos,1)); }
			
			//cerr << md[p] << "|";
			p_old = p;
			p = md.find_first_of("ACGT^",p_old+1);
		}
	}
	num = md.substr(p_old+1);				
	//cerr << num << endl;
}

// Fasta_Read
// By default, it finds all microsatellites that are at least 8bp long (total length), 
// where the repeat sequence is between 1bp and 4bp, and is repeated at least 3 times.

// $ (echo ">A"; echo "CAAAAAAAAAAAAAAAAAAAAAAAT") > a.fa
// $ ./quickscan -f a.fa
// Processing sequences in a.fa...
// A len=25
// 2 24 23 A 22+1 cAAAAAAAAAAAAAAAAAAAAAAAt 0
// finished in 8.9e-05s
//////////////////////////////////////////////////////////////////////////
/*
void findTandems(const string & seq, const string & tag)
{
	
	const int OFFSET_TABLE_SIZE = 100;

	int MIN_REPORT_LEN = 8;
	int MIN_REPORT_UNITS = 3;
	int MIN_UNIT_LEN = 1;
	int MAX_UNIT_LEN = 4;
	int FLANK = 10;
	
	cout << ">" << tag << " len=" << seq.length() << endl;

	//if (MUTATE_fp != NULL)
	//{
	//	fprintf(MUTATE_fp, ">%s MUTATE_PROB=%0.03f MUTATE_MAX_LEN=%d MUTATE_PROB_EXTEND=%0.03f MUTATE_SEED=%d\n", 
	//	tag.c_str(), MUTATE_PROB, MUTATE_MAX_LEN, MUTATE_PROB_EXTEND, MUTATE_SEED);
	//}
	
	int lastprinted = 0;

	int offsets[OFFSET_TABLE_SIZE][OFFSET_TABLE_SIZE];

	// initialize the offsets
	for (int merlen = 1; merlen <= MAX_UNIT_LEN; merlen++)
	{
		for (int phase = 0; phase < merlen; phase++)
		{
			offsets[merlen][phase] = phase;
		}
	}

	// now scan the sequence, considering mers starting at position i
	for (int i = 0; i < seq.length(); i++)
	{
    
		// consider all possible merlens from 1 to max
		for (int merlen = 1; merlen <= MAX_UNIT_LEN; merlen++)
		{
			int phase = i % merlen;
			int offset = offsets[merlen][phase];

			// compare [i..i+merlen) to [offset..offset+merlen)
			int j = 0;
			while ((j < merlen) && 
				(i+j < seq.length()) && 
					(seq[i+j] == seq[offset+j])) 
						{ j++; }

			// is the end of the tandem?
			if (j != merlen || (i+j+1 == seq.length()))
			{
				// am i the leftmost version of this tandem?
				if (seq[offset-1] != seq[offset+merlen-1])
				{
					// is it long enough to report?
					if (((i-offset)/merlen >= MIN_REPORT_UNITS) && (i - offset >= MIN_REPORT_LEN))
					{
						// is it primitive?
						int ml = 1;

						while (ml < merlen)
						{
							int units = (i-offset+j) / ml;

							int allmatch = 1;
							for (int index = 1; allmatch && (index < units); index++)
							{
								// compare the bases of the current unit to those of unit0
								for (int m = 0; m < ml; m++)
								{
									if (seq[offset+m] != seq[offset+index*ml+m])
									{
										allmatch = 0;
										break;
									}
								}
							}

							if (!allmatch) { ml++; }
							else           { break; }
						}

						// everything checks, now report it

						if (ml == merlen)
						{
							// start end length
							cout << offset+1 << "\t" << i+j << "\t" << i+j-offset << "\t";

							// tandem seq
							for (int z = 0; z < merlen; z++) { cout << seq[offset+z]; }

							// complete units + remainder
							cout << "\t" << (i - offset) / merlen << "+" << j << "\t";

							// left flank - tandem - right flank
							for (int z = offset-FLANK; z < offset;    z++) { if (z >= 0) { cout << (char) tolower(seq[z]); } }
							for (int z = offset;       z < i+j;       z++) { cout << seq[z]; }
							for (int z = i+j;          z < i+j+FLANK; z++) { if (z < seq.length()) { cout << (char) tolower(seq[z]); } }

							cout << "\t" << phase;


							if (MUTATE_fp != NULL)
							{
								float r = ((float) rand()) / ((float) RAND_MAX);

								if (r <= MUTATE_PROB)
								{
									int mutate_len = 1 + (rand() % MUTATE_MAX_LEN);
									float extend = ((float) rand()) / ((float) RAND_MAX);
									char mutate = '-';
									if (extend <= MUTATE_PROB_EXTEND) { mutate = '+'; }

									cout << "\t||\t" << mutate << "\t" << mutate_len << "\t";

									// left flank
									for (int z = offset-FLANK; z < offset;  z++) { if (z >= 0) { cout << (char) tolower(seq[z]); } }

									// mutated ms
									if (mutate == '-')
									{
										for (int z = offset; z < i+j-mutate_len; z++) { cout << seq[z]; }
									}
									else
									{
										for (int z = offset; z < i+j; z++) { cout << seq[z]; }
										for (int xxx = 0; xxx < mutate_len; xxx++) { cout << seq[offset + ((j+xxx)%merlen)]; }
									}

									// right flank
									for (int z = i+j; z < i+j+FLANK; z++) { if (z < seq.length()) { cout << (char) tolower(seq[z]); } }


									// genome sequence up to ms
									for (int ppp = lastprinted; ppp < offset; ppp++) { fprintf(MUTATE_fp, "%c", seq[ppp]); }

									// mutated ms
									if (mutate == '-')
									{
										for (int z = offset; z < i+j-mutate_len; z++) { fprintf(MUTATE_fp, "%c", seq[z]); }
									}
									else
									{
										for (int z = offset; z < i+j; z++) { fprintf(MUTATE_fp, "%c",  seq[z]); }
										for (int xxx = 0; xxx < mutate_len; xxx++) { fprintf(MUTATE_fp, "%c", seq[offset + ((j+xxx)%merlen)]); }
									}

									lastprinted = i+j;
								}
							}

							cout << endl;
						}
					}
				}

				offsets[merlen][phase] = i;
			}
		}
	}
}
*/