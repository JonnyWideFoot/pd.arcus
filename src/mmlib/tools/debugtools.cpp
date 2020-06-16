#include "global.h"

#if _DEBUG

#include "maths/maths.h"
#include "workspace/workspace.h"
#include "DebugTools.h"

void DebugTools::printDebugAxis( int &atomNumberStart, int &molNumber )
{
	printDebugAtomLine( atomNumberStart++, 'Q', molNumber, 0, 0, 0 );
	printDebugAtomLine( atomNumberStart++, 'Q', molNumber, 1, 0, 0 );
	printDebugAtomLine( atomNumberStart++, 'Q', molNumber, 0, 1, 0 );
	printDebugAtomLine( atomNumberStart++, 'Q', molNumber, 0, 0, 1 );
	printf("\n");
	molNumber++;
}

void DebugTools::printDebugOrigin( int &atomNumberStart, int &molNumber )
{
	printDebugAtomLine( atomNumberStart++, 'Q', molNumber, 0, 0, 0 );
	printf("\n");
	molNumber++;
}

#endif

