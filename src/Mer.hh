#ifndef MER_HH
#define MER_HH 1

/******************************************************************
** Mer.hh
**
** Class for storing the representation of a canonical Mer: 
** contiguos sub-sequence of bases from one input read
**
**  Authors: Giuseppe Narzisi & Michael C. Schatz
**    Date: December 11, 2013
**
*******************************************************************/

#include <string>
#include "util.hh"


using namespace std;

typedef char Ori_t;
const Ori_t F = 'F';
const Ori_t R = 'R';

// Mer_t
//////////////////////////////////////////////////////////////////////////

typedef string Mer_t;

// CanonicalMer_t
//////////////////////////////////////////////////////////////////////////

class CanonicalMer_t
{
public:

	Mer_t mer_m;
	Ori_t ori_m;

	CanonicalMer_t() {}

	CanonicalMer_t(Mer_t mer) 
		{ set(mer); }

	void set(Mer_t mer)
	{
		Mer_t rmer = rc2(mer);

		if (mer < rmer)
		{
			mer_m = mer;
			ori_m = F;
		}
		else
		{
			mer_m = rmer;
			ori_m = R;
		}
	}

	ostream & print(ostream & out) const
	{
		//out << mer_m << ":" << ((ori_m == F) ? 'F' : 'R');
		string str = ((ori_m == F) ? ":F" : ":R");
		out << mer_m << str;
		return out;
	}

	// CanonicalMer_t print
	/////////////////////////////////////////////////////////////////

	friend ostream& operator<<(std::ostream& o, const CanonicalMer_t & mer)
	{
		return mer.print(o);
	}

	static Mer_t rc(const Mer_t & mer)
	{
		Mer_t retval;

		for (int i = mer.length()-1; i >= 0; i--)
		{
			retval.push_back(rrc(mer[i]));
		}

		return retval;
	}
	
	
	static Mer_t rc2(const Mer_t & mer) {
	    unsigned char c;
 	   
	   	int i = 0;
		int j = mer.length()-1;
	   	Mer_t retval = mer;
		
		while(i<=j) {
 
			c = rrc(retval[i]);
			retval[i]= rrc(retval[j]);
			retval[j] = c;
			i++; j--;
	    }
		
		return retval;
	}
	
};


#endif
