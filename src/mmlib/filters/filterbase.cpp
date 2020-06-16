#include "global.h"
#include "filterbase.h"

FilterBase::FilterBase()
{
}

FilterBase::~FilterBase()
{
}

void FilterBase::setInternalName()
{
	name = "FilterBase";
}

MolFilterBase::MolFilterBase() 
	: FilterBase(), 
	m_Mol(NULL) 
{
}

MolFilterBase::MolFilterBase( const MoleculeBase& _mol ) 
	: FilterBase(), 
	m_Mol(&_mol) 
{
}

void MolFilterBase::setInternalName()
{
	name = "MolFilterBase";
}

FilterContainer::FilterContainer()
{
}

void FilterContainer::setInternalName()
{
	name = "FilterContainer";
}

bool FilterContainer::passes()
{
	for( size_t i = 0; i < size(); i++ )
	{
		if( !element(i).passes() ) 
			return false;
	}
	return true;
}

std::string FilterContainer::reason()
{
	StringBuilder sb;
	bool foundBad = false;
	for( size_t i = 0; i < size(); i++ )
	{
		if( !element(i).passes() ) 
		{
			if( !foundBad ) 
			{
				sb.appendFormat("Filter Collection `%s' failed:\n")(name);
				foundBad = true;
			}
			sb.appendFormat("  (%d/%d) `%s' failed: %s\n")
				(i+1)
				(size())
				(element(i).name.c_str())
				(element(i).reason().c_str());						
		}
	}
	if( !foundBad ) sb.appendFormat("All %d Internal filters passed\n")(size());
	return sb.toString();
}

