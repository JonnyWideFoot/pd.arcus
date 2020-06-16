#ifndef __NUMSASA_H
#define __NUMSASA_H

#include "workspace/workspace.fwd.h"






//-------------------------------------------------
//
/// \brief  Implements calculation of numerical SASAs
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API NumSASA
{
public:
	NumSASA();
	NumSASA( const std::string& _DatFileName );

	// read out data file
	void readDat(const std::string& _DatFileName);
	void setTo( const WorkSpace& _WSpace );

	// Info
	void detail() const;
	void info() const;

	/// Calculate the SASA
	void calc();
	
private:
	void assertData();

	//std::vector<dvector> m_Sphere;
};

#endif


