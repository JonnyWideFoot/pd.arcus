#include "global.h"
#include "tools/io.h"
#include "library/libpath.h"

LibraryPathStore* LibraryPathStore::getSingleton()
{
	static LibraryPathStore inst;
	return &inst;
}

LibraryPathStore::LibraryPathStore()
{
	parseEnvironmentVariable();
}

// Parses the environment variable PD_PARAM_PATH into search paths to look for
// library files. The order of searching is: local  directory, lib path,
void LibraryPathStore::parseEnvironmentVariable()
{
	m_LibraryPath.clear();
	m_LibraryPath.push_back(".");

	char *c_env = getenv("PD_PARAM_PATH");

	// if the return value of getenv is NULL the variable is undefined.
	if(c_env == NULL)
	{
		printf("Library Paths:  PD_PARAM_PATH undefined\n");
		return;
	}

	std::string env_PD_PARAM_PATH = "";
	env_PD_PARAM_PATH = c_env;
	std::vector <std::string> temp_m_LibraryPath = chopstr(env_PD_PARAM_PATH, ";:, \t\12\15" );

	for(size_t i=0;i<temp_m_LibraryPath.size();i++)
	{
		m_LibraryPath.push_back( temp_m_LibraryPath[i] );
	}

	printf("Library Paths: \n");
	for(size_t i=0;i<m_LibraryPath.size();i++)
	{
		printf("  %s/\n",m_LibraryPath[i].c_str());
	}
}

std::string LibraryPathStore::findFullFilename(const std::string &filename)
{
	for(int i=0;i<m_LibraryPath.size();i++)
	{
		std::string composite;
		composite = m_LibraryPath[i] + "/" + filename;
		if(IO::fileExists(composite))
		{
			return composite;
		};
	}
	return filename; // return filename evwen if it fails such that error messages are fine
}

