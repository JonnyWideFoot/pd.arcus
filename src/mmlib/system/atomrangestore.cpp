#include "global.h"
#include "system/molecule.h"
#include "atomrangestore.h"

// -----------------
//  AtomRangeStore
// -----------------

AtomRangeStore::AtomRangeStore() : m_Molecule(NULL)
{
}

AtomRangeStore::AtomRangeStore( MoleculeBase &_Molecule, const PickAtomRange& _Range ) : m_Molecule(NULL)
{
	setTarget(_Molecule,_Range);
}

void AtomRangeStore::setTarget( MoleculeBase &_Molecule, const PickAtomRange& _Range )
{
	m_Molecule = &_Molecule;
	m_Range = _Range;
	m_PosCache.resize(m_Range.getNAtoms()); // will refill the vector with dvector() (default constructor yields 0,0,0)
}

void AtomRangeStore::store()
{
	size_t count = m_PosCache.size();
	size_t start = m_Range.getStartAtomIndex();

	ASSERT( m_Molecule != NULL, CodeException, "AtomRangeStore::store() called when setTarget() has not been called!");
	D_ASSERT( (count == (m_Range.getEndAtomIndex()-start+1)), CodeException,
		"AtomRangeStore::store() internal cache is the wrong size!");

	for( size_t i = 0; i < count; i++ )
	{
		m_PosCache[i].setTo(m_Molecule->atomxyz(start+i));
	}
}

void AtomRangeStore::revert()
{
	size_t count = m_PosCache.size();
	size_t start = m_Range.getStartAtomIndex();

	ASSERT( m_Molecule != NULL, CodeException, "AtomRangeStore::revert() called when setTarget() has not been called!");
	D_ASSERT( (count == (m_Range.getEndAtomIndex()-start+1)), CodeException,
		"AtomRangeStore::revert() internal cache is the wrong size!");

	for( size_t i = 0; i < count; i++ )
	{
		m_Molecule->atomxyz(start+i).setTo(m_PosCache[i]);
	}
}

