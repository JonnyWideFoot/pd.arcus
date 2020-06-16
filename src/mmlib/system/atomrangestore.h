#ifndef _ATOM_RANGE_STORE
#define _ATOM_RANGE_STORE

#include "pickers/basicpickers.h"

class MoleculeBase;

//-------------------------------------------------
/// \brief  A specialisation of a PosStore for an atomic range as opposed to an
/// atomic picker. There may be efficiency to be had by doing this...
/// The class was written before PosStore , which is capable of the same thing with similar syntax,
/// but requires pointer lookups to store, it may therefore be slower. (unbenchmarked)
/// We should discuss wether to keep or remove this class.
/// \author Jon Rea 
class AtomRangeStore
{
public:
	AtomRangeStore();
	AtomRangeStore( MoleculeBase &_Molecule, const PickAtomRange& _Range ) ;
	void setTarget( MoleculeBase &_Molecule, const PickAtomRange& _Range );
	void store();
	void revert();

protected:
	std::vector<Maths::dvector> m_PosCache;
	PickAtomRange m_Range;
	MoleculeBase* m_Molecule;
};

#endif



