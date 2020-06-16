#ifndef __ARCUS_BASE_H
#define __ARCUS_BASE_H

#include "protocols/protocolbase.h" // Provides the 'ProtocolBase' base class
#include "manipulators/conformerbuilder.h" // Provides the 'ConfBuild_RandomSingle' member
#include "filters/filterbase.h" // Provides 'FilterContainer'
#include "workspace/workspace.fwd.h"

// Forward Declarations
class PickResidueRange;

namespace Physics
{
	class FF_BreakableBonded;
	class IBondBreakable;
}

namespace Library
{
	class AngleSet;
}

namespace Protocol
{
	//-------------------------------------------------
	//
	/// \brief A rapid loop stitching method onto a protein core
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	///
	/// \author Jon Rea Copyright 2006/7
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API ArcusBase : public ProtocolBase
	{
	public: 
		/// Default constructor
		ArcusBase( Physics::Forcefield& _ffs, Physics::Forcefield& _ff, const Library::AngleSet& _as );

		bool ExtraReport; ///< Outputs more information. e.g. The initial structure is appended to the tra file.
		bool FlagRebuilt; ///< Should we deem the output from this class as rebuilt? If true, we set Particle::isRebuildRequired() to false following build.
		bool UtiliseIdealised; ///< Use the idealised state for execution, or what we started with?
		bool MinimiseIdealised; ///< Minimise the bond angle torsion and improper terms for the segment prior to execution?
		bool SetRefStatePostInitialRebuild; ///< When the internal minor-atom builders are invoked, should that be the reference state (required for H's)

		// Filtering
		std::string SegFanFilename; ///< If this filename is set, the SegDistFan will be enabled

		void buidAutoDetect(); ///< Detect which regions to build from the isRebuildRequired() Particle flag
		void buildAddExplicit( const PickResidueRange& Range ); // explicitly define a region to build
		void clear(); ///< Clean up any region definitions and set to a fully reinitialsed state.

		enum FinalStateType
		{
			LastAcc, 
			LowestEpot
		} FinalState;

		/// prints a little block of parameter information
		virtual void info() const; 		

		virtual int run();
		virtual int runcore() = 0; 

		const Library::AngleSet& getAngleSet() const { return *m_AngleSet; }

	protected:

		virtual int initialise(); ///< Called in the intitialisation phase of run_core()

		size_t Step;
		Physics::Forcefield* ffs;

		bool m_SplitBuilder; ///< Some derived-classes require that there are 2 builders per region, if dual-branch.

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       

		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

		std::vector<SegmentDef> m_Regions; ///< This holds the definitions for the sections which must be built.
		std::vector<SegmentDef> m_ConfBuilderSeg; ///< Required to prevent the SegmentDef going out of scope, as the conformer builder depends on it
		std::vector<Manipulator::ConfBuild_RandomSingle> m_ConfBuilder;
		std::vector<FilterContainer> m_Filters;

		void setPickedRegions( const WorkSpace& wspace, PickAtomRanges& picker, std::vector<SegmentDef> _regions ) const;
		void setPickedRegions( const WorkSpace& wspace, PickAtomRanges& picker, SegmentDef& _region ) const;
		void setPickedRegions_Core( const WorkSpace& wspace, PickAtomRanges& picker, SegmentDef& _region ) const;		

	private:
		bool verifyNoRegionClash(); ///< make sure none of the build regions overlap

		virtual int addFFBreaks(); ///< Break the relevent bonds, returning the number broken
		int addFFBreaks( std::vector<Physics::IBondBreakable*>& _bondedComponents );

		bool initialBuild(); ///< Create the regions in an extended conformation.
		virtual void configLevel1Filters( size_t i );

		void callFlagRebuiltIf(); ///< If FlagRebuilt is flagged, set isRebuildRequired() to false for all built atoms.
		void setFlagRebuilt( bool setBuilt ); ///< Set the rebuild state for the atoms in the current residue ranges

		std::vector<Physics::IBondBreakable*> m_BondedFFS;
		std::vector<Physics::IBondBreakable*> m_BondedFF; ///< This is assigned during runcore(). A pointer to the underlying breakable forcefield is required for ARCUS to function.
		const Library::AngleSet* m_AngleSet;

		bool m_AutoRegions;
	};
}

#endif

