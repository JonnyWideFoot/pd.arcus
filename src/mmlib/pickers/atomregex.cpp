#include "global.h"
#include "tools/stringbuilder.h"
#include "system/fundamentals.h"
#include "atomregex.h"

// ------------------------------------
// Begin Class PickAtomRegEx
// ------------------------------------

PickAtomRegEx::PickAtomRegEx( const std::string& _pattern ) 
{ 
	setPattern(_pattern); 
}

PickAtomRegEx* PickAtomRegEx::clone() const
{
	return new PickAtomRegEx(*this);
}

void PickAtomRegEx::setPattern( const std::string& _pattern )
{
	m_AtomNames.clear();
	m_Pattern = _pattern;
	StringBuilder sb( _pattern );
	size_t index;
	while( (index = sb.LastOf(PickAtomRegExDelimiter)) != SIZE_T_FAIL )
	{
		m_AtomNames.push_back(sb.toString(index,sb.size()-index));
		sb.erase(index,sb.size()-index);
	}
	m_AtomNames.push_back(sb.toString());
}

bool PickAtomRegEx::matches( const Particle& particle ) const
{
	for( size_t i = 0; i < m_AtomNames.size(); i++ )
	{
		if( 0 == particle.pdbname.compare( m_AtomNames[i] ) ) 
		{
			return true;
		}
	}
	return false;
}

