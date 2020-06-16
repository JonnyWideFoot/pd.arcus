#ifndef __ATOM_REGEX_H
#define __ATOM_REGEX_H

#include <string>
#include <vector>

#include "pickers/pickbase.h" // PickBase base class

const char PickAtomRegExDelimiter = '-';

//-------------------------------------------------
/// \brief  A simple RegEx atom picker
/// \details Currently this is the simplest of the simple REGEX for atoms, only individual names can be
/// expressed, separated by the 'PickAtomRegExDelimiter'.
/// This class should be extended to allow much more complex atom and parent residue expressions
/// but at the moment I don't have the time to code it...
/// \author Jon Rea 
/// \todo State of Development: Primordial Soup
/// \bug More than likely
class PickAtomRegEx : public PickBase
{
public:
	// Public constructor logic
	PickAtomRegEx();
	PickAtomRegEx( const std::string& _pattern );
	virtual PickAtomRegEx* clone() const;

	// Public function calls
	void setPattern( const std::string& _pattern ); ///< set the atom selection patten
	virtual bool matches( const Particle& particle ) const; ///< Returns true if the passed reference matches the internal pattern
protected:
	std::string m_Pattern;
	std::vector<std::string> m_AtomNames;
};

#endif

