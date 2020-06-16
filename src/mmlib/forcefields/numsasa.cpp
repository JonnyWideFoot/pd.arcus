#include "global.h"
#include <fstream>
#include "tools/stringtool.h"
#include "pickers/pickbase.h"
#include "workspace/workspace.h"
#include "numsasa.h"

NumSASA::NumSASA()
{
}

NumSASA::NumSASA( const std::string& _DatFileName )
{
	readDat(_DatFileName);
}

void NumSASA::assertData()
{
	//ASSERT( m_WSpace != NULL, ProcedureException, "The internal workspace pointer is null; setTo() has not been called");
	//ASSERT( m_Sphere.size() > 0, ProcedureException, "The Numeric SASA has not been initialised with its sphere data");
	//ASSERT( m_AtomIndexes.size() > 0, ProcedureException, "There is no sasa data; has calc been called? or does the  
}

void NumSASA::readDat(const std::string& _DatFileName)
{
	//m_Sphere.clear();
}

void NumSASA::setTo( const WorkSpace& _WSpace )
{
	//m_Picker = PickAllAtoms();
}

void NumSASA::detail() const
{
}

void NumSASA::info() const
{
}

void NumSASA::calc()
{
}


