#ifndef __FASTA_H
#define __FASTA_H

#include "sequence/sequence.h" // provides the base class






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// Index this collection with a simple string obtained after 
/// the '>' on the first line of each sequence definition
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API Fasta : public Sequence::BioSequenceCollection<std::string> 
{
public:
	Fasta();
	virtual ~Fasta();

	virtual void parseLine( const std::string &_line );
};

#endif

