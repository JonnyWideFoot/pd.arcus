#include "global.h"
#include "filters/filterbase.h"
#include "manipulators/conformerbuilderbase.h"
#include "workspace/segdef.h"
#include "conformer.h"

using namespace Manipulator;

ConformerFileWriter::ConformerFileWriter() :
	m_File(NULL),
	m_Begun(false),
	m_Source(NULL),
	m_Filter(NULL),
	m_WriteCount(0)
{
}

ConformerFileWriter::~ConformerFileWriter()
{
	CloseHandle();
}

void ConformerFileWriter::begin( const std::string& _Filename, ConfBuilderBase& _Donor )
{
	if( m_Begun ) throw ProcedureException("ConformerFileWriter is already open. Call end() before opening again!");
	m_Begun = true;

	// Set our conformer source
	m_Source = &_Donor;

	// Attempt to open the file
	m_File = fopen(_Filename.c_str(),"w");
	if(m_File == NULL) throw IOException("Could not open file for ConformerFileWriter!");

	// Print sequence header
	const Sequence::BioSequence& seq = m_Source->getSegDef().getSequence();
	std::string seqS = seq.printToStringSingle();
	fprintf( m_File, "%s\n", seqS.c_str() );
}

bool ConformerFileWriter::push_back()
{
	if( !m_Begun ) throw ProcedureException("ConformerFileWriter is not open. Call begin() before push_back()!");
	m_Begun = false;

	if( m_Filter != NULL && !m_Filter->passes() ) 
	{
		// The filter does not accept this structure
		return false;
	}
	else
	{
		// We either have no filter enabled, or the structure passes.
		size_t loopLength = m_Source->getSegDef().getNRes();
		const Conformer conf = m_Source->getConformer();
		for( size_t i = 0; i < loopLength; i++ )
		{
			// 48 is used to convert the small positive number to an ASCII character
			fprintf( m_File, "%c", (char)(conf[i] + 48) );
		}
		fprintf( m_File, "\n");
		m_WriteCount++;
		return true;
	}
}

void ConformerFileWriter::end()
{
	if( !m_Begun ) throw ProcedureException("ConformerFileWriter is not open. Call begin() before closing again!");
	m_Begun = false;

	// Re-initialsise to null
	m_File = NULL;
	m_Source = NULL;
	m_Filter = NULL;
	m_WriteCount = 0;
}

void ConformerFileWriter::CloseHandle()
{
	fclose(m_File);
}

