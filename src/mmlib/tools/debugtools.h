#if _DEBUG

#ifndef DEBUG_TOOLS_H
#define DEBUG_TOOLS_H

#include "maths/maths.fwd.h"

class PD_API Particle;

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API DebugTools
{
public:
	static void DebugTools::printDebugAxis( int &atomNumberStart, int &molNumber );
	static void DebugTools::printDebugOrigin( int &atomNumberStart, int &molNumber );
	static void DebugTools::printDebugAtomLine( Particle *atom, int index );

	template <class T>
	static void printDebugAtomLine( int AtomNumber, char element, int MolNumber, Maths::Tvector<T> *v )
	{
		printDebugAtomLine( AtomNumber, element, MolNumber, v->x, v->y, v->z );
	}

	template <class T>
	static void printDebugAtomLine( int AtomNumber, char element[4], int MolNumber, Maths::Tvector<T> *v )
	{
		printDebugAtomLine( AtomNumber, element, MolNumber, v->x, v->y, v->z );
	}

	template <class T>
	static void printDebugAtomLine(int AtomNumber, char element[4], int MolNumber, T x, T y, T z )
	{
		printf("ATOM %5d %c%c%c%c ALA %3d %8.3lf%8.3lf%8.3lf%\n",
			AtomNumber,
			element[0],
			element[1],
			element[2],
			element[3],
			MolNumber,
			x, y, z);
	}

	template <class T>
	static void printDebugAtomLine(int AtomNumber, char element, int MolNumber, T x, T y, T z )
	{
		printf("ATOM  %5d  %c   MOL   %3d    %8.3f%8.3f%8.3f%\n",
			AtomNumber,
			element,
			MolNumber,
			x, y, z);
	}


private:
	DebugTools()
	{
	}
};

#endif // DEBUG_TOOLS_H

#endif // _DEBUG
