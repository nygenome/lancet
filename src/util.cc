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


// return file name without extension
StringType GetBaseFilename(const char *filename)
{
    StringType fName(filename);
    size_t pos = fName.rfind(".");
    if(pos == StringType::npos)  //No extension.
        return fName;

    if(pos == 0)    //. is at the front. Not an extension.
        return fName;

    return fName.substr(0, pos);
}

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
string dtos(double d) // convert double to string
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

/*
------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
A             A           Adenine
C             C           Cytosine
G             G           Guanine
T             T           Thymine
U             U           Uracil
M           A or C
R           A or G
W           A or T
S           C or G
Y           C or T
K           G or T
V         A or C or G
H         A or C or T
D         A or G or T
B         C or G or T
X or N  G or A or T or C
N       G or A or T or C
. or -       gap
*/

// isAmbiguos 
// return true if the base is an ambiguos IUPAC code
//////////////////////////////////////////////////////////////////////////

bool isAmbiguos(char b)
{
	if (b == 'M' || b == 'm' ||
		b == 'R' || b == 'r' ||
		b == 'W' || b == 'w' ||
		b == 'S' || b == 's' ||
		b == 'Y' || b == 'y' ||
		b == 'K' || b == 'k' ||
		b == 'V' || b == 'v' ||
		b == 'H' || b == 'h' ||
		b == 'D' || b == 'd' ||
		b == 'B' || b == 'b' ||
		b == 'X' || b == 'x')
	{ 
		return true;
	}

	return false;
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

	for (int i = str.length()-1; i >= 0; --i)
	{
		retval.push_back(rrc(str[i]));
	}

	return retval;
}

// reverse
// reverse the content of a string (in place)
//////////////////////////////////////////////////////////////
void reverse(string & str) 
{
	string retval = str;
	
	int i=0;
	int j=str.size()-1;
	while(i<j){
		char tmp = str[i];
		str[i] = str[j];
		str[j] = tmp; 
		++i;--j;
	}	
}

// isNseq
// returns true if the input sequence contains only Ns
//////////////////////////////////////////////////////////////////////////
bool isNseq(const string & seq)
{
	bool result = true;

	int end = seq.length();
	for (int i = 0; i < end; ++i)
	{

		if ( (seq[i] != 'N') || (seq[i] != 'n') ) {
			result = false;
			break;
		}
	}
	return result;	
}

// HammingDistance
// returns the hamming distance between two strings or -1 if strings have different length
//////////////////////////////////////////////////////////////////////////
int HammingDistance(const string & s1, const string & s2)
{
	int dist = 0;

	if (s1.length() != s2.length() ) { return -1; }
	
	for (unsigned int i = 0; i < s1.length(); ++i) {
		if(s1[i] != s2[i]) { dist++; }
	}
	
	return dist;
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
	for (int offset = 0; offset < end; ++offset)
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
	for (int offset = 0; offset < end; ++offset)
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
				++count;
				if(count>max) { flag = false; break; }
			}
			++j;
		}
		//if(flag && j==p.size()) { return true; }
		if(flag && j==L) { return true; }
		++i;
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
void parseMD(string & md, map<int,int> & M, int start, string & qual, int min_qv) {
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
			++pos;
			assert( (pos-start) >= 0 );
			if(qual[pos-start] >= min_qv) {
				//cerr << qual[pos-start] << ">=?" << min_qv << endl;
				mit = M.find(pos);
				if (mit != M.end()) { ++((*mit).second); }
				else { M.insert(std::pair<int,int>(pos,1)); }
			}
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
bool findTandems(const string & seq, const string & tag, int max_unit_len, int min_report_units, int min_report_len, int dist_from_str, int pos, int & len, std::string & motif)
{
	bool ans = false;
	//FILE * MUTATE_fp = NULL;
	const int OFFSET_TABLE_SIZE = 100;

	unsigned int MIN_REPORT_LEN = min_report_len;
	unsigned int MIN_REPORT_UNITS = min_report_units;
	//unsigned int MIN_UNIT_LEN = 1;
	unsigned int MAX_UNIT_LEN = max_unit_len;
	unsigned int FLANK = 10;
	int delta = dist_from_str;
	
	//cerr << ">" << tag << " len=" << seq.length() << endl;

	//if (MUTATE_fp != NULL)
	//{
	//	fprintf(MUTATE_fp, ">%s MUTATE_PROB=%0.03f MUTATE_MAX_LEN=%d MUTATE_PROB_EXTEND=%0.03f MUTATE_SEED=%d\n", 
	//	tag.c_str(), MUTATE_PROB, MUTATE_MAX_LEN, MUTATE_PROB_EXTEND, MUTATE_SEED);
	//}
	
	//int lastprinted = 0;

	int offsets[OFFSET_TABLE_SIZE][OFFSET_TABLE_SIZE];

	// initialize the offsets
	for (unsigned int merlen = 1; merlen <= MAX_UNIT_LEN; ++merlen)
	{
		for (unsigned int phase = 0; phase < merlen; ++phase)
		{
			offsets[merlen][phase] = phase;
		}
	}

	// now scan the sequence, considering mers starting at position i
	for (unsigned int i = 0; i < seq.length(); ++i)
	{
    
		// consider all possible merlens from 1 to max
		for (unsigned int merlen = 1; merlen <= MAX_UNIT_LEN; ++merlen)
		{
			int phase = i % merlen;
			int offset = offsets[merlen][phase];

			// compare [i..i+merlen) to [offset..offset+merlen)
			unsigned int j = 0;
			while ((j < merlen) && 
				(i+j < seq.length()) && 
					(seq[i+j] == seq[offset+j])) 
						{ ++j; }

			// is the end of the tandem?
			if (j != merlen || (i+j+1 == seq.length()))
			{	
				assert(offset-1 < (int)seq.length());
				assert(offset+merlen-1 < seq.length());
											
				// am i the leftmost version of this tandem?
				if (seq[offset-1] != seq[offset+merlen-1])
				{					
					// is it long enough to report?
					if (((i-offset)/merlen >= MIN_REPORT_UNITS) && (i - offset >= MIN_REPORT_LEN))
					{
						// is it primitive?
						unsigned int ml = 1;

						while (ml < merlen)
						{
							unsigned int units = (i-offset+j) / ml;

							int allmatch = 1;
							for (unsigned int index = 1; allmatch && (index < units); ++index)
							{
								// compare the bases of the current unit to those of unit0
								for (unsigned int m = 0; m < ml; ++m)
								{
									if (seq[offset+m] != seq[offset+index*ml+m])
									{
										allmatch = 0;
										break;
									}
								}
							}

							if (!allmatch) { ++ml; }
							else           { break; }
						}

						// everything checks, now report it

						if (ml == merlen)
						{
							// start end length
							//cerr << offset+1 << "\t" << i+j << "\t" << i+j-offset << "\t";
							
							int start = offset;
							int end = i+j;
							int L = i+j-offset;
							
							if ( (pos >= (start-delta)) && (pos <= (end+delta)) ) { 
								ans = true; 
								// store STR motif and size
								len = L; 
								for (unsigned int z = 0; z < merlen; ++z) { 
									motif += seq[offset+z];
								}
							}

							// tandem seq
							/*
							for (unsigned int z = 0; z < merlen; z++) { 
								cerr << seq[offset+z];
							}
							*/
							// complete units + remainder
							//cerr << "\t" << (i - offset) / merlen << "+" << j << "\t";

							// left flank - tandem - right flank
							for (unsigned int z = offset-FLANK; z < (unsigned int)offset;    ++z) { if (z >= 0) { /*cerr << (char) tolower(seq[z]);*/ } }
							for (unsigned int z = offset;       z < i+j;       ++z) { /*cerr << seq[z];*/ }
							for (unsigned int z = i+j;          z < i+j+FLANK; ++z) { if (z < seq.length()) { /*cerr << (char) tolower(seq[z]);*/ } }

							//cerr << "\t" << phase;

							/*
							if (MUTATE_fp != NULL)
							{
								float r = ((float) rand()) / ((float) RAND_MAX);

								if (r <= MUTATE_PROB)
								{
									int mutate_len = 1 + (rand() % MUTATE_MAX_LEN);
									float extend = ((float) rand()) / ((float) RAND_MAX);
									char mutate = '-';
									if (extend <= MUTATE_PROB_EXTEND) { mutate = '+'; }

									cerr << "\t||\t" << mutate << "\t" << mutate_len << "\t";

									// left flank
									for (int z = offset-FLANK; z < offset;  z++) { if (z >= 0) { cerr << (char) tolower(seq[z]); } }

									// mutated ms
									if (mutate == '-')
									{
										for (unsigned int z = offset; z < i+j-mutate_len; z++) { cerr << seq[z]; }
									}
									else
									{
										for (unsigned int z = offset; z < i+j; z++) { cerr << seq[z]; }
										for (unsigned int xxx = 0; xxx < mutate_len; xxx++) { cerr << seq[offset + ((j+xxx)%merlen)]; }
									}

									// right flank
									for (unsigned int z = i+j; z < i+j+FLANK; z++) { if (z < seq.length()) { cerr << (char) tolower(seq[z]); } }


									// genome sequence up to ms
									for (unsigned int ppp = lastprinted; ppp < offset; ppp++) { fprintf(MUTATE_fp, "%c", seq[ppp]); }

									// mutated ms
									if (mutate == '-')
									{
										for (unsigned int z = offset; z < i+j-mutate_len; z++) { fprintf(MUTATE_fp, "%c", seq[z]); }
									}
									else
									{
										for (unsigned int z = offset; z < i+j; z++) { fprintf(MUTATE_fp, "%c",  seq[z]); }
										for (unsigned int xxx = 0; xxx < mutate_len; xxx++) { fprintf(MUTATE_fp, "%c", seq[offset + ((j+xxx)%merlen)]); }
									}

									lastprinted = i+j;
								}
							}
							*/
							//cerr << endl;
						}
					}
				}

				offsets[merlen][phase] = i;
			}
		}
	}
	return ans;
}
