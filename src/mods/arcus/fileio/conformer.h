#ifndef __CONFORMER_IO_H
#define __CONFORMER_IO_H

#include "typedefs.h"

// Forward declarations
namespace Manipulator
{
	class ConfBuilderBase;
}
class FilterBase;

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class ConformerFileWriter
{
private:
	ConformerFileWriter( const ConformerFileWriter& _copy ); ///< Copy operators are private. As a file handle is stored, copying of this class is not allowed!
	ConformerFileWriter& operator=( const ConformerFileWriter& _copy );

public:
	ConformerFileWriter();
	~ConformerFileWriter();

	void begin( const std::string& _Filename, Manipulator::ConfBuilderBase& _Donor ); ///< Inform this class of its source and open a file handle.
	bool push_back(); ///< Add the next conformer from the source given to Begin(), but only if the current filter passes true. Returns the push result.
	void end(); ///< Complete the write process and close the file handle
	void SetFilter( FilterBase& _filter ) { m_Filter = &_filter; }
	void RemoveFilter() { m_Filter = NULL; }
	conformer_count_type getWriteCount() const { return m_WriteCount; }

private:
	// Member Helper Functions
	void CloseHandle();

	// Member Data
	FILE* m_File;
	bool m_Begun;
	Manipulator::ConfBuilderBase* m_Source;
	FilterBase* m_Filter;
	conformer_count_type m_WriteCount;
};

#endif

