#ifndef __REVISION
#define __REVISION

#include <string>

//-------------------------------------------------
/// \brief Holds PDs Revision information
/// \author Mike Tyka & Jon Rea
class Revision
{
private:
	Revision(); // You cant create an instance of this class
public:
	static std::string LongRevString();
	static std::string ShortRevString();
	static std::string LibVersion();
	static std::string LibName();
	static std::string FullStamp();
	static std::string CopyRight();
private:
	static void doInit();
	static bool init;
	static std::string longRev;
	static std::string shortRev;
	static std::string fullStamp;
};

#endif

