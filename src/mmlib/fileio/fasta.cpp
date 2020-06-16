#include "global.h"
#include "sequence/sequence.h"
#include "fasta.h"

using namespace Sequence;

Fasta::Fasta() : BioSequenceCollection<std::string>()
{
	THROW(NotImplementedException,"");
}

Fasta::~Fasta()
{
}

void Fasta::parseLine( const std::string &_line )
{
	THROW(NotImplementedException,"");
}

