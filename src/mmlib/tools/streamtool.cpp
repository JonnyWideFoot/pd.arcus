#include "global.h"

#include "streamtool.h"

#include <iomanip> // Needed for some stream manipulation functionality e.g. setw()

using namespace std;

void setFloatFmt( std::ostream &_file, int _width, int _precision )
{
	_file.setf(ios::showpoint); // ensure the decimal point is set
	_file.setf(ios::right,ios::adjustfield); // ensure all is right aligned to the available width
	_file.setf(ios::fixed,ios::floatfield); // ensure that the floating point representation
	_file.width(_width); // set the print width
	_file.precision(_precision); // set the number of decimal palaces
}

