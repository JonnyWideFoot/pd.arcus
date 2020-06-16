#include "global.h"
#include "maths/maths.h"
#include "system/molecule.h"
#include "system/rebuilder.h"
#include "pickers/atomregex.h"
#include "segdef.h"

// -------------
// SegmentDefBase
// -------------

SegmentDefBase::SegmentDefBase() :
	PickResidueRange(),
	m_WSpace(NULL),
	m_Init(false),
	m_NativeKnown(false),
	m_NativeBBKnown(false),
	m_BreakType(SegBreakStart)
{
}

SegmentDefBase::SegmentDefBase( WorkSpace& _wspace, int _startindex, int _endindex, SegBreakType _break ) :
	PickResidueRange(),
	m_WSpace(NULL),
	m_Init(false),
	m_NativeKnown(false),
	m_NativeBBKnown(false),
	m_BreakType(SegBreakStart) // Set in the call to init()
{
	init(_wspace,_startindex,_endindex,_break);
}

size_t SegmentDefBase::getBreakResIndex() const
{
	switch( m_BreakType )
	{
	case SegBreakCentre:
		{
			return getCentreResIndex()+1;
		}
	case SegBreakEnd:
		{
			return getEndResIndex()+1;
		}
	case SegBreakStart:
		{
			return getStartResIndex();
		}
	default:
		THROW(CodeException,"Unknown SegBreakType encountered!");
	}
}

void SegmentDefBase::init( WorkSpace& _wspace, int _startindex, int _endindex, SegBreakType _break )
{
	// Derived classes depend upon immutability of this class to allow performance releated memory caching - we therefore only allow one initialisation.
	ASSERT( !m_Init, CodeException, "SegmentRange init() has already been called!");
	m_Init = true; // flag this - we can only ever do this once to ensure the immutabiluty of this class instance

	m_BreakType = _break;

	// Set the range within the base clas - this is validated 
	PickResidueRange::setRange(_wspace,_startindex,_endindex);
	
	// Set the internal pointer
	m_WSpace = &_wspace;

	// setup native calls
	m_NativeKnown = true;
	m_NativeBBKnown = true;
	for( size_t i = getStartResIndex(); i <= getEndResIndex(); i++ )
	{
		int start = m_WSpace->res[i].ifirst;
		int end = m_WSpace->res[i].ilast;
		for( int j = start; j <= end; j++ )
		{
			if( !m_WSpace->atom[j].isKnownStructure() )
			{
				if( !m_WSpace->atom[j].isHydrogen() )
				{
					// We only care if its a heavy atom, the hydrogens are almost exclusively
					// rebuilt for every experimenatal structure.
					m_NativeKnown = false;

					if( m_WSpace->atom[j].isBackbone() )
					{
						// If we dont know the non-hydrogen-backbone atoms, we dont know the native state
						m_NativeBBKnown = false;
					}
				}				
				if( !m_NativeKnown && !m_NativeBBKnown )
				{
					goto DoneNativeFind; // Break the outer loop, both are false...
				}
			}
		}
	}

DoneNativeFind:
	// setup our internal sequence
	m_Sequence.append(m_WSpace->getSequence(),getStartResIndex(),getNRes());
}

Particle& SegmentDefBase::operator[]( size_t _index )
{
	D_ASSERT( _index < m_NAtoms, OutOfRangeException, "Atom out of range of Segment!");
	return m_WSpace->atom[_index + m_StartAtomIndex];
}

const Particle& SegmentDefBase::operator[]( size_t _index ) const
{
	D_ASSERT( _index < m_NAtoms, OutOfRangeException, "Atom out of range of Segment!");
	return m_WSpace->atom[_index + m_StartAtomIndex];
}

Particle& SegmentDefBase::at( size_t _index )
{
	D_ASSERT( _index < m_NAtoms, OutOfRangeException, "Atom out of range of Segment!");
	return m_WSpace->atom[_index + m_StartAtomIndex];
}

const Particle& SegmentDefBase::at( size_t _index ) const
{
	D_ASSERT( _index < m_NAtoms, OutOfRangeException, "Atom out of range of Segment!");
	return m_WSpace->atom[_index + m_StartAtomIndex];
}

bool SegmentDefBase::isNativeKnown() const
{
	return m_NativeKnown;
}

bool SegmentDefBase::isNativeBackboneKnown() const
{
	return m_NativeBBKnown;
}

WorkSpace& SegmentDefBase::getWorkSpace() const
{
	return *m_WSpace;
}

const Sequence::BioSequence& SegmentDefBase::getSequence() const
{
	return m_Sequence;
}

void SegmentDefBase::info() const
{
	std::string seq = getSequence().printToStringSingle();
	Printf("Segment Info:\nLength: %d\nSequence: %s\nBreak-Class: %d\n")
		(getNRes())(seq.c_str())((int)m_BreakType);
}

SegStore::SegStore()
{
	m_DataContained = false;
}

SegStore::SegStore( SegmentDefBase& _mol, const PickBase& _picker ) 
{
	setPicking( _mol, _picker );
}

void SegStore::setPicking( SegmentDefBase& _seg, const PickBase& _picker )
{	
	m_WSpace = &_seg.getWorkSpace();
	m_Start = _seg.getStartAtomIndex();
	m_Length = _seg.getNAtoms();
	lookupAtoms(_picker);
	m_DataContained = false;
}


SegmentDef::SegmentDef() :
	SegmentDefBase()
{
}

SegmentDef::SegmentDef( WorkSpace& _mol, int _startindex, int _length, SegBreakType _break ) :
	SegmentDefBase(_mol,_startindex,_length,_break)
{
	initCore( _mol, _startindex, _length, _break );
}

void SegmentDef::init( WorkSpace& _mol, int _startindex, int _length, SegBreakType _break )
{
	// Ensure we initialise the base class too. 
	// The base class will throw an excepton if init has already occured.
	SegmentDefBase::init(_mol,_startindex,_length,_break);
	// Now the core init for this class.
	initCore( _mol, _startindex, _length, _break );
}

void SegmentDef::initCore( WorkSpace& _mol, int _startindex, int _length, SegBreakType _break )
{
	// Proxy to the underlying base class
	SegmentDefBase& segDef = *this;

	//m_BBNoCAIndex.setPicking(segDef,PickAtomRegEx("C-O-N"));

	m_CAStore.setPicking(segDef,PickAtomRegEx("CA"));
	m_CAStore.store();
	m_BBStore.setPicking(segDef,PickAtomRegEx("CA-C-O-N"));
	m_BBStore.store();
	m_CBStore.setPicking(segDef,PickAtomRegEx("CB"));
	m_CBStore.store();
	m_HeavyStore.setPicking(segDef,PickHeavyAtoms());
	m_HeavyStore.store();

	m_Original.setTarget(_mol,segDef);
	m_Original.store();

	m_Cache.setTarget(_mol,segDef);	
	m_Cache.store();
}

//const PosPointer& SegmentDef::getBackboneNoCAPositions() const
//{
//	return m_BBNoCAIndex;
//}

