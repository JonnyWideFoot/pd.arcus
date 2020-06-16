#ifndef _LOOP_DEFINITION_H
#define _LOOP_DEFINITION_H

#include "workspace/pospointer.h" // Defines the position caching members
#include "fileio/tra.h" // Defines the TraBlock base classes
#include "sequence/sequence.h" // Defines the BioSequence class member
#include "pickers/basicpickers.h" // Defines the PickResidueRange base class
#include "system/atomrangestore.h" // Provides a class member
#include "maths/fastrandom.h"
#include "workspace/workspace.fwd.h"

class MoleculeBase;
class Particle;

enum SegBreakType
{
	SegBreakEnd,
	SegBreakStart,	
	SegBreakCentre
};

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
/// Class defining a continuous section of a given molecule and the implied atomic ranges.
///
/// \details DETAILED USER'S DESCRIPTION
/// Accessors for the sequence and wether the native structure is defined are also provided.
/// This class caches atomic ranges, and therefore is incompatible with MoleculeBase, which 
/// can change size after allocation. We therefore must deal only with the Immutable 
/// workspace class.
///
/// \author Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class SegmentDefBase : public PickResidueRange
{
public:
	SegmentDefBase();
	SegmentDefBase( WorkSpace& _mol, int _startindex, int _endindex, SegBreakType _break );
	void init( WorkSpace& _mol, int _startindex, int _endindex, SegBreakType _break );
	bool isInitialised() const { return m_Init; }
	
	// Manipulation Orientation
	inline SegBreakType getBreakType() const { return m_BreakType; }
	size_t getBreakResIndex() const; //< Returns the residue whos Omega angle is broken by torsional changes in this Seg

	// Particles and positions
	Particle& at( size_t _index ); ///< Access the particle at the zero based _index from the start of this segment
	const Particle& at( size_t _index ) const; ///< Access the const particle at the zero based _index from the start of this segment
#ifndef SWIG
	Particle& operator[]( size_t _index ); ///< Access the particle at the zero based _index from the start of this segment
	const Particle& operator[]( size_t _index ) const; ///< Access the const particle at the zero based _index from the start of this segment
#endif

	// Native?
	bool isNativeKnown() const; ///< checks to see isKnownStructure is set for all heavy atoms.
	bool isNativeBackboneKnown() const; ///< checks to see isKnownStructure is set for all backbone atoms.
	
	// Parent
	WorkSpace& getWorkSpace() const; ///< Returns a reference to the parent molecule

	// Sequence
	const Sequence::BioSequence& getSequence() const; ///< Returns the sequence of this segment

	void info() const;

private:
	WorkSpace* m_WSpace;
	Sequence::BioSequence m_Sequence;

	bool m_NativeKnown; ///< Is the native heavy atom structure known for this segment? isNative terms are claculated on initialisation and cached here.
	bool m_NativeBBKnown; ///< Is the native backbone structure known for this segment? isNative terms are claculated on initialisation and cached here.

private:
	SegBreakType m_BreakType; ///< Should the conformer be manipulated from the forward or reverse anchor?
	bool m_Init; ///< // Record initialisation privatly - we can only ever do this once to ensure the immutabiluty of each SegmentDefBase class instance
};


class PD_API SegStore : public PosStore
{
public:
	SegStore();
	SegStore( SegmentDefBase& _mol, const PickBase& _picker );

	/// set a new position pointer
	void setPicking( SegmentDefBase& _seg, const PickBase& _picker );
};


//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class SegmentDef : public SegmentDefBase, public IO::BTF_Block_Comment_User, public IO::BTF_Block_Vector_User
{
public:
	SegmentDef(); ///< Waits for init() to be called at a later time
	SegmentDef( WorkSpace& _mol, int _startindex, int _endindex, SegBreakType _break ); ///< Internally calls init() following construction
	void init( WorkSpace& _mol, int _startindex, int _endindex, SegBreakType _break ); ///< Allows optional post-constructor initialisation

	//const PosPointer& getCAPositions() const;
	//const PosPointer& getBackbonePositions() const;
	//const PosPointer& getBackboneNoCAPositions() const;
	//const PosPointer& getCBPositions() const;
	
	double cCRMS_CA()    { return m_CAStore.calcCRMS(); }    ///< CA atom cRMS of the current structure against the reference state
	double cCRMS_BB()    { return m_BBStore.calcCRMS(); }    ///< Backbone atom cRMS of the current structure against the reference state
	//double cCRMS_AA()    { return m_AAStore.calcCRMS(); }    ///< All atom cRMS of the current structure against the reference state
	double cCRMS_Heavy() { return m_HeavyStore.calcCRMS(); } ///< All heavy atom cRMS of the current structure against the reference state

	inline void storeCache()  { m_Cache.store(); } ///< Store the current state of the segment using the internal cache.
	inline void revertCache() { m_Cache.revert(); } ///< Revert to the segment to the current cached state
	inline void revertOriginal() { m_Original.revert(); } ///< Revert the segment to the structure defined when this class was initialised.

protected:

	// Store initial positions from RMS calculations
	SegStore m_CAStore;
	SegStore m_BBStore;
	SegStore m_CBStore;
	SegStore m_HeavyStore;

	AtomRangeStore m_Original;

	/// A dynamic store with public accessors that can cache a position set at any given time.
	/// Example usage would include post-idealisation of a segment during a building process.
	AtomRangeStore m_Cache;

private:
	void initCore( WorkSpace& _mol, int _startindex, int _endindex, SegBreakType _break );
};


//-------------------------------------------------
//
/// \brief  
/// SegmentDefUser derviced classes function on an internal 'SegmentDef', which encapsulated the
/// features of a given seg within the member particlesystem.
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
class PD_API SegmentDefUser
{
public:
	SegmentDefUser( SegmentDef& _segdef ) :
		m_SegLength( _segdef.getNRes() ),
		m_SegStartResIndex( _segdef.getStartResIndex() ),
		m_SegCentreResIndex( _segdef.getCentreResIndex() )
	{
		ASSERT( _segdef.isInitialised(), CodeException, "SegmentDef is not initialised!");
		m_SegDef = &_segdef; // assign the protected member.
		m_WSpace = &_segdef.getWorkSpace(); // retrieve the pointer to the loops parent 
	}
	
	inline SegmentDef& getSegDef() const { return *m_SegDef; } ///< Public accessor to the underlying full segment definition
	inline WorkSpace& getWSpace() const { return *m_WSpace; } ///< Public accessor to the underlying workspace

	inline size_t segLength() const { return m_SegLength; } ///< Public read-only accessor for the length of this defined segment
	inline size_t segStartIndex() const { return m_SegStartResIndex; } ///< Public read-only accessor for the start index of this defined segment
	inline size_t segCentreIndex() const { return m_SegCentreResIndex; } ///< Public read-only accessor for the start index of this defined segment

private:
	SegmentDef* m_SegDef; ///< Private member pointer to the underlying full segment definition
	WorkSpace* m_WSpace; ///< Private member pointer to the underlying workspace
	size_t m_SegLength; ///< Private cached index values. This represents the length of the represented segment.
	size_t m_SegStartResIndex; ///< Private cached index values. This represents the start index of the represented segment.
	size_t m_SegCentreResIndex; ///< Private cached index values. This represents the centre index of the represented segment.
};






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
class AnchorBase
{
};






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
class PeptideAnchor : public AnchorBase
{
};

#endif


