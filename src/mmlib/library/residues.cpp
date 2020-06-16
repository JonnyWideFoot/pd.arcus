#include "global.h"

#include "residues.h"

namespace Library
{
	// --------------------------------------------------------------------------------
	// Enumeration of standard residue types and their grouped definitions
	// --------------------------------------------------------------------------------

	StandardResidues getResidue(char _QureyResidue)
	{
		// Protein
		if(_QureyResidue=='A') return A;
		else if(_QureyResidue=='C') return C;
		else if(_QureyResidue=='D') return D;
		else if(_QureyResidue=='E') return E;
		else if(_QureyResidue=='F') return F;
		else if(_QureyResidue=='G') return G;
		else if(_QureyResidue=='H') return H;
		else if(_QureyResidue=='I') return I;
		else if(_QureyResidue=='K') return K;
		else if(_QureyResidue=='L') return L;
		else if(_QureyResidue=='M') return M;
		else if(_QureyResidue=='N') return N;
		else if(_QureyResidue=='P') return P;
		else if(_QureyResidue=='p') return p;
		else if(_QureyResidue=='Q') return Q;
		else if(_QureyResidue=='R') return R;
		else if(_QureyResidue=='S') return S;
		else if(_QureyResidue=='T') return T;
		else if(_QureyResidue=='V') return V;
		else if(_QureyResidue=='W') return W;
		else if(_QureyResidue=='Y') return Y;

		// Nucleotide
		else if(_QureyResidue=='a') return a;
		else if(_QureyResidue=='t') return t;
		else if(_QureyResidue=='g') return g;
		else if(_QureyResidue=='c') return c;
		else if(_QureyResidue=='u') return u;

		// Oh well...
		else return Unknown;
	}

	bool HasResidueProperty(StandardResidueProperties _Includes, StandardResidues _QureyResidue)
	{
		return ((unsigned int)_Includes & (unsigned int)_QureyResidue) > 0;
	}

	bool HasResidueProperty(StandardResidueProperties _Includes, char _QureyResidue)
	{
		return ((unsigned int)_Includes & (unsigned int)getResidue(_QureyResidue)) > 0;
	}

	bool HasSameResidueProperty( StandardResidues _QureyResidue1, StandardResidues _QureyResidue2 )
	{
		if( ((_QureyResidue1 & ShortHydrophobic) > 0) && ((_QureyResidue2 & ShortHydrophobic) > 0) ) return true;
		if( ((_QureyResidue1 & BulkyAromatic) > 0) && ((_QureyResidue2 & BulkyAromatic) > 0) ) return true;
		if( ((_QureyResidue1 & Polar) > 0) && ((_QureyResidue2 & Polar) > 0) ) return true;
		if( ((_QureyResidue1 & Pyrimidine) > 0) && ((_QureyResidue2 & Pyrimidine) > 0) ) return true;
		if( ((_QureyResidue1 & Purine) > 0) && ((_QureyResidue2 & Purine) > 0) ) return true;
		if( ((_QureyResidue1 & Positive) > 0) && ((_QureyResidue2 & Positive) > 0) ) return true;
		if( ((_QureyResidue1 & Negative) > 0) && ((_QureyResidue2 & Negative) > 0) ) return true;
		return false;
	}

	bool HasSimilarResidueProperty( StandardResidues _QureyResidue1, StandardResidues _QureyResidue2 )
	{
		THROW(NotImplementedException,"");
	}

	bool HasSameResidueProperty( char _QureyResidue1, char _QureyResidue2 )
	{
		StandardResidues r1 = getResidue(_QureyResidue1);
		if( r1 == Unknown ) return false;
		StandardResidues r2 = getResidue(_QureyResidue2);
		if( r2 == Unknown ) return false;
		return HasSameResidueProperty(r1,r2);
	}

	bool HasSimilarResidueProperty( char _QureyResidue1, char _QureyResidue2 )
	{
		StandardResidues r1 = getResidue(_QureyResidue1);
		if( r1 == Unknown ) return false;
		StandardResidues r2 = getResidue(_QureyResidue2);
		if( r2 == Unknown ) return false;
		return HasSameResidueProperty(r1,r2);
	}
}

