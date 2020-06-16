#ifndef _PICK_FROM_FILE
#define _PICK_FROM_FILE

#include <string>
#include <vector>
#include "sequence/sequence.h"
#include "pickers/basicpickers.h"




//-------------------------------------------------
//
/// \brief  This class is needed by PickResiduesFromFile.
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
class ResidueList
{
public:
	ResidueList();

#ifndef SWIG
	inline const Sequence::BioSource& operator[](int i) { return m_ResList[i]; };
#endif

	void setChainID(const char& _chainID); ///< set the ChainID for this residue
	char getChainID() const; ///< get the ChainID for this residue

	void pushback(const Sequence::BioSource& _resInfo); ///< add an instance of a BioSource to the std::vector m_ResInfo
	int getResListSize() const; ///< returns m_ResList.size()
	int getBiologicalIndex(int _ResListIndex) const;
	char getBiologicalIcode(int _ResListIndex) const;
	std::string getResName(int _ResListIndex) const;

protected:
	char m_ChainID;                             ///< store of ChainID.
	std::vector<Sequence::BioSource> m_ResList; ///< store of Biological Index, Biological Icode and Residue Name for each residue.
};






//-------------------------------------------------
//
/// \brief  This class defines a list of residues supplied in an input file.
///
/// \details 
/// Use a seperate instance of the class (and a seperate input file) for each molecule.
/// The input file should contain a list of residues, each on a new line containing the
/// information (in this order, seperated by whitespace):
/// chain ID, residue number, residue 3 letter code (or *) and residue Icode (optional).
/// # can be used to comment out a line.
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PickResiduesFromFile: public PickResidueBase
{
public:
	PickResiduesFromFile(const MoleculeBase& _molBase); ///< Default constructor if user needs the picker to exist but doesn't want to actually pick any residues.
	PickResiduesFromFile(const MoleculeBase& _molBase, const std::string& _filename); ///< Constructor. Each residue listed in the file should be on a new line.

	virtual PickResiduesFromFile* clone() const;

	virtual bool matches( const Particle& particle ) const;
	virtual bool matches( const Residue& _res ) const;

protected:	
	std::vector<ResidueList> m_ResListByChain;
	std::vector<bool> m_ResIsPicked; ///< size can be set to the number of residues in the workspace.

	/// Read in the list of residues from a file.
	/// Each residue listed in the file should be on a new line
	/// and each entry should contain (in this order):
	/// a chain ID, a residue number, a residue 3 letter code and a residue Icode (optional).
	void loadList(const std::string& _filename);

	void compareAll(const MoleculeBase& _molBase); ///< initialisation function for the matches() function. Compare m_ResListByChain to the information in the workspace.
};

#endif

