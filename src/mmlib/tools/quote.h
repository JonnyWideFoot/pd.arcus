#ifndef __QUOTE_H
#define __QUOTE_H

#include <string>
#include <vector>






//-------------------------------------------------
//
/// \brief A Quote Engine :-) 
///
/// \details No MM library would be complete without a quote engine ;-)
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API Quote
{
public:
	static Quote* getSingleton(); /// 'Singleton Pattern' Implementation
	Quote();
	void IncludeDB( std::string _filename );
	void printQuote();
private:
	std::vector<std::string> m_Quotes;
};

#endif

