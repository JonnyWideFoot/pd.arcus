// Pre-compiled Header
#include "global.h"

#include "tools/stringtool.h"

// Self Header Include Should Be Last
#include "example.h"

namespace Physics
{
	void exampleStandAloneFunction()
	{
		std::string test = "teST";		
		makeupper(test);
		printf("String: %s \n", test.c_str());
	}
} 

