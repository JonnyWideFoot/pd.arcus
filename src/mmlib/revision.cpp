#include "global.h"
#include "tools/stringbuilder.h"
#include "revision.h"

bool Revision::init = false;
std::string Revision::longRev = "";
std::string Revision::shortRev = "";
std::string Revision::fullStamp = "";

void Revision::doInit()
{
	// Auto-make the subversion version string.
	StringBuilder build( "$Revision: 1184 $" ); // *NOTE* this string is automatically edited by subversion when revision.cpp is altered
	build.Trim("$"); // remove $ signs
	build.erase(0,9); // remove 'Revision:'
	build.Trim(); // remove remaining whitespace

	shortRev = build.toString();
	longRev = "SVN Revision: " + shortRev;

	build.clear();

	build.append( LibName() );
	build.append( ' ' );
	build.append( LibVersion() );
	build.append( ' ' );
	build.append( longRev );

	fullStamp = build.toString();

	init = true;
}

std::string Revision::FullStamp()
{
	if( !init ) doInit();
	return fullStamp;
}

std::string Revision::CopyRight()
{
	return "(c) J.Rea & M.Tyka 2003-2007";
}

std::string Revision::LibName()
{
	return "MMLib";
}

std::string Revision::LibVersion()
{
	return "v0.8b";
}

std::string Revision::LongRevString()
{
	if( !init ) doInit();
	return longRev;
}

std::string Revision::ShortRevString()
{
	if( !init ) doInit();
	return shortRev;
}

