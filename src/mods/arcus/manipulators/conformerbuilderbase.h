#ifndef __CONFORMER_BUILDER_BASE_H
#define __CONFORMER_BUILDER_BASE_H

// system includes
#include <vector>

// mmlib includes
#include "verbosity.h"
#include "typedefs.h" // Required for the use of the 'conformer_type' typedef
#include "library/backbonetorsions.h" // Required for member 'BackboneTorsionSubSet'
#include "workspace/segdef.h" // Required for 'SegmentDefUser' base class
#include "segtorsions.h" // Required for 'SegCoreBBOnlyTC' and 'SegWholeBackboneTC' base classes

// Forward delclarations
#include "maths/fastrandom.fwd.h"

namespace Library
{
	class AngleSet;
}

namespace Manipulator
{
	typedef std::vector<conformer_type> Conformer;

	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// \note
	///  The main 4 + A single residue specific sidechain atom. (not-implemented)
	/// I now think that a better way to implement this would be to calculate the centroid on the fly from the
	/// backbone coords.
	/// Advantages:
	/// 1) Which atom is used for the 'centroid' is best left to the FF i.e. CB for most H2 for GLY - delegating
	///    to the protocol means that the conf builder would need to read the FFPS.
	/// 2) Exactly where to place the centroid, e.g. a given average distance along the CA-CB vector?,
	///    is also better left to the FF during its calc.
	/// 3) Any given centroid FF can also be implemented as a FilterBase and added to a normal conformer
	///    enumeration.
	/// 4) The concept of the centroid is abstracted away from the protocol. The conformer builder should by
	///    definition really only produce either a backbone or all atom structure conforming to a given angleset.
	/// 5) Doing so scales linearly with conformer length. Whereas Moving extra CB/H's around scales much less
	///    well, although there is a tradeoff between the two with varying length.
	///    Overall I dont think this is the primary consideration.
	///
	/// \author Mike Tyka & Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class ConfBuilderBase :
		public virtual SegmentDefUser,
		public SegCoreBBOnlyTC,
		public SegWholeBackboneTC
	{

	public:
		ConfBuilderBase(
			const Library::AngleSet& _angleSet, ///< The angleset that will be used by the conformer builder
			SegmentDef &_segdef ///< The defintion for the section of the parent molecule that will be considered
			);

		// ----------
		// build mode
		// ----------
		static const int ConformerBuildMode_Count = 2;
		enum ConformerBuildMode
		{
			Backbone = 0, ///< Just the main 4 backbone atoms (N, CA, C, O)
			AllAtom = 1 ///< All the available atoms, including dummies
		};

		void MinimiseIdealised(); ///< Should the idealised structure be 'energy' minimised on bonds angles and torsions only; Defaults to true.
		Verbosity::Type OutputLevel; ///< STFU ?

		// Do we want the builer to add random noise around a given angle-pair? 
		// How much by? This value sets the 1.0 s.d. point of a normally distributed delta value
		// moving away from the actual conformer angle.
		double NoiseSigma;

		// ------------
		//  Build Mode
		// ------------
		void setBuildMode( ConformerBuildMode mode ); ///< use this to change wether the builder is using rapid backbone or allatom mode. Stage1 generation should occur under backbone only mode.
		inline ConformerBuildMode getBuildMode(){ return m_BuildMode; }

		// ------------------
		// Internal Conformer
		// ------------------
		inline const Conformer& getConformer() { return m_CurrentConformer; }
		inline conformer_count_type getIndex() { return m_CurrentIndex; }
		bool isCurrentConformerNative();

		// ------------------------
		//  Structure Manipulation
		// ------------------------
		inline void revertIdealised() { m_Idealised.revert(); applyWholeConformer(); }
		inline void revertClosestToNative() { m_ClosestToNative.revert(); }
		virtual void reset() = 0; ///< provides mode specific reinitialisation
		virtual bool next() = 0; ///< increment to the next valid conformer and manipulate the particle system. Returns true if a conformer is available, false if not.

		// --------------------
		// public user feedback
		// --------------------
		virtual void info() const; ///< sequence, anglecounts, native info if defined
		inline double getBestPossibleARMS() const { return m_BestARms; }
		void printCurrentConformer(); ///< user feedback of the current conformer state
		void printCurrentInternalConformer(); ///< used only in debuging
		std::string sprintCurrentConformer(); ///< print the current conformer to a string
		std::string sprintCurrentInternalConformer(); ///< used only in debuging
		void printAngleSubset() const;

		// ----------------------------------------------
		// cRMS calls.
		// We have both cA cRMS and allBackboneAtoms cRMS
		// ----------------------------------------------
		inline double getCurrentConformationCRMS_CA() { return getSegDef().cCRMS_CA(); } // CRMS: Alpha Carbon
		inline double getCurrentConformationCRMS_BB() { return getSegDef().cCRMS_BB(); } // CRMS: All non-hydrogen backbone atoms
		inline void printCurrentConformationCRMS_CA() { printf("%lf", getCurrentConformationCRMS_CA() ); }
		inline void printCurrentConformationCRMS_BB() { printf("%lf", getCurrentConformationCRMS_BB() ); }

		inline const Library::AngleSet& getAngleSet(){ return *m_AngleSet; }

	protected:
		// whats the best we can do with the available angleset?
		bool assessNearestToNativeConformer(); ///< Find the native conformer, log it, and append to the BristolTrajectoryFormat file as a sanity check...
		double m_BestARms; ///< the best aRMS structure for the current set of conformers within the current angleset

		/// A random numbr generator for internal use..
		Maths::FastRandom* m_Rand;

		// angleset pointers
		const Library::AngleSet* m_AngleSet; ///< external pointer/memory cleanup
		std::vector<Library::BackboneTorsionSubSet> m_BackboneTorsionSets; ///< Array: one per loop residue
		void addEntireAnglesetLibrary(); ///< For each loop residue import all the angleset states from the m_AngleSet into the m_BackboneTorsionSets members
		void cloneBackboneTorsionSubSet( ConfBuilderBase& _Builder  ); ///< Set the current internal angleset to that of another builder

		// Conformer data :
		//
		// ---------- !! UBER SPECIAL NOTES, TAKE NOTICE !! -----------
		//
		// m_CurrentConformer holds the ResidueAngle array index of the* *REDUCED** angle subset, **NOT** the original
		// index in the main '.angleset' library
		// m_NativeConformer however holds the array index in the '.angleset' library!
		// The conformers are stored as 'conformer_type' arrays (unsigned short, defined in 'conformertype.h')
		ConformerBuildMode m_BuildMode;
		Conformer m_CurrentConformer;
		Conformer m_NativeConformer;
		std::string m_NativeNumChar; ///< the native conformer characters i.e. '1' not 1, '2' not 2
		std::string m_NativeBinChar; ///< the native conformer-bin characters
		conformer_count_type m_CurrentIndex; ///< the current conformer in the sequence defined by the enumeration method.

		// Structure Manipulation
		void applyWholeConformer(); ///< apply the whole current conformer based on the current build mode
		void performRotation( int conformerIndex ); ///< perform a given loop structure roation based on the current build mode

	private:
		// Private helper functions
		void verifyAngleDefs(); ///< Make sure the sequence of the segment that we have is available in the internal angleset
		void performIdealisation(); ///< Internal call to idealise the internal segment definiton - i.e. make all the bonds and angles 'forcefield ideal' from the forcefield position definitions

		// Private Data
		AtomRangeStore m_Idealised; ///< A segstore to cache the idealised coordinates
		AtomRangeStore m_ClosestToNative; ///< A segstore to cache the coordinates of the closest conformer to the native state in terms of aRMS

		bool m_IdealisedMinimised;
	};

	enum DescriptorType
	{
		AllResidues,
		NativeBin,
		NativeClosest
	};


	//-------------------------------------------------
	//
	/// \brief  
	/// A base class for all true conformer builders. Takes a descriptor defining the underlying angle subset,
	/// which defaults on initialisation to the whole library, but can be changed after construction.
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
	class ConfBuilderBase_Descriptor: public ConfBuilderBase
	{
	public:
		ConfBuilderBase_Descriptor(
			const Library::AngleSet &angleSet, ///< Our chosen angleset
			SegmentDef &segdef ///< The protein segment on which random perturbation should occur.
			);

		void setDescriptor( const std::string& confDescriptor );
		void setDescriptor( DescriptorType _mode );

		conformer_count_type getPermutations() const { return m_PossibleConformers; }

	protected:
		void infoDescriptor() const;		

	private:
		void validateDescriptor();
		void calcTotalAngles(); ///< multiply up the number in the angle subsets
		void initDescriptor();
		void obtainBBSet();

		conformer_count_type m_PossibleConformers; ///< potentially a BIG number
		std::string m_Descriptor;
	};


	//-------------------------------------------------
	//
	/// \brief  A conformer system that is enumerated using a random perturbation process
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
	class ConfBuilderBase_Random: public ConfBuilderBase_Descriptor
	{
	public:
		ConfBuilderBase_Random(
			const Library::AngleSet &angleSet, ///< Our chosen angleset
			SegmentDef &segdef ///< The protein segment on which random perturbation should occur.,
			);

		void setRandCount( conformer_count_type randomCount, bool warnOverMax = true );
		conformer_count_type getRandCount() const { return m_RandConformers; }
		virtual void reset();

		bool PropensityWeighting; ///< Should random changes be weighted by propensity?

		int getPreviousChangeIndex() const { return m_ChangeIndex; }

	protected:
		int m_ChangeIndex;
		void RandomiseWholeConformer();
		void RandomiseSingleConformerPos(); ///< A random change at a random index
		void RandomiseSingleConformerPos( size_t _changeIndex ); ///< A randon change at a specific index
		conformer_count_type m_RandConformers;
	};


	//-------------------------------------------------
	//
	/// \brief   A conformer system that is enumerated by sequential additions to internal conformer positions
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
	class ConfBuilderBase_Enum: public ConfBuilderBase_Descriptor
	{
	public:
		enum ConformerEnumMode
		{
			NumericAscent, ///< Increment the conformer in a column (like a normal number) - more simple and rapid
			SequentialAscent ///< Increment the conformer in a row - less simple and rapid, but allows better ordering by propensity if the conformer list is pre-sorted (it is)
		};

		ConfBuilderBase_Enum(
			Library::AngleSet &angleSet, ///< Our chosen angleset
			SegmentDef &segdef ///< The protein segment on which perturbation should occur.
			);

		void changeEnumerationMode( ConformerEnumMode _Mode ); ///< Change the way that the current conformer is enumerated

		virtual void reset();
		virtual bool next();

	protected:
		ConformerEnumMode m_EnumerationMode;

		//NumericAscent
		void incrementDescriptorAndapply( int conformerIndex ); ///< called by the next() function, recursive function to increment and then apply the new conformer (low geometry-change calculation waste, nice :-D)

		//SequentialAscent
		void initStates();
		bool nextState(); ///< generate the m_State_roots, before 1st call to nextSubState()
		void initSubState();
		bool nextSubState(); ///< enumerate valid posibilities generated GenerateStates
		bool isValid(); ///< Used to prevent generation of two posibilities for '1123' (i.e. 1a1b23 and 1b1a23) and to check the conformer indexes generated are within range.
		void printStateDefinition();
		void applySequentialAscentConformer();

		Conformer m_State_MaxIndx;
		conformer_type m_State_MaxMax;
		Conformer m_State_Root;
		Conformer m_State_Current;

		int m_StateResultantIndex;
		std::vector<int> m_State_RepeatCounter;
		std::vector<int> m_SortTemplate_Current;
		std::vector<int> m_SortTemplate_Root;

		std::vector<bool> m_SubStateTaken;
		std::vector<int> m_SubState_HadLast; //int needed to store an integer index, including a -1 'NULL' value
		int m_SubStatePosition;

	private:
		void ClearSequentialModeMemory(); ///< Free up some memory when we are not using it
	};


	//-------------------------------------------------
	//
	/// \brief  A conformer system that is enumerated by sequential selection from a given conformer buffer
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
	class ConfBuilderBase_Buffer: public ConfBuilderBase
	{
	public:
		ConfBuilderBase_Buffer(
			const Library::AngleSet &angleSet, ///< Our chosen angleset
			SegmentDef &segdef ///< The protein segment on which perturbation should occur.
			);

		bool push_back( const Conformer& _Conformer ); ///< add an explicit conformer to the buffer (used by the FromFile() conformerbuilder mode)

		virtual void reset();
		virtual bool next();

		inline void clear() { m_Conformers.clear(); } /// clear all the conformers out of the current buffer
		inline size_t capacity() const { return m_InternalCapacity; }
		inline size_t size() const { return m_Conformers.size(); }

		void setInfinateCapacity() { m_InternalCapacity = 0; } ///< Set the internal capacity value to 0, meaning infinate.
		void setMaxCapacity( size_t maxBufferCapacity );

		void printBuffer(); ///< print the entire buffer to the standard output. Use with **caution** if very large ;-)

	protected:
		bool rotateFromBuffer( conformer_count_type conformerIndex ); ///< set the ConformersCurrent state to the conformer in the buffer
		bool validateConformer( const Conformer& _Conformer ); ///< We require validation before appending to the buffer.
		std::vector< Conformer > m_Conformers; ///< Dynamic array to hold the conformers.
		size_t m_InternalCapacity; ///< Limit conformer import into the above vector, 0 to disable.
	};
}

#endif

