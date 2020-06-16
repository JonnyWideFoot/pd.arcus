#ifndef __ARCUS_REFINER_H
#define __ARCUS_REFINER_H

#include "fileio/trablocks.h"
#include "protocols/minimise.h"
#include "forcefields/restraintbase.h"
#include "forcefields/restraint_positional.h"
#include "manipulators/rotamer_scwrl.h"
#include "arcus_base.h"

namespace Library
{
	class RotamerLibrary;
}

namespace Protocol
{
	/// \brief This is what we will be refining...
	struct LoopSet
	{
		PickResidueRange range;
		std::vector<PosStore> posCache;
		std::vector< std::pair<double,PosStore*> > posCachePointer;
		std::vector<size_t> clusterReps; ///< Indexes in the **SORTED** posCachePointer, not the posCache!!!
	};

	/// \brief Defines a class that can perform loop refinement
	/// \details ARCUS takes a reference to a refiner, so it can be overridden
	/// All refiners derive from ArcusRefineBase.
	class PD_API ArcusRefineBase : public IO::BTF_Block_Comment_User
	{
	public:
		ArcusRefineBase() : AppendBestToTra(false), m_WSpace(NULL), m_FF(NULL), OutputLevel(Verbosity::Normal) {}
		void initialise( WorkSpace& _wspace, Physics::Forcefield& _ff );
		virtual void refine( std::vector<LoopSet>& loopSet ) = 0; ///< Refine all clusterReps
		bool AppendBestToTra;
		WorkSpace& getWSpace() const { return *m_WSpace; }

		Verbosity::Type OutputLevel;

	protected:
		WorkSpace* m_WSpace;
		Physics::Forcefield* m_FF;
	};

	/// \brief Simpler algorithms refining each loop in isolation
	class PD_API ArcusRefineIndividual : public ArcusRefineBase
	{
	public:
		ArcusRefineIndividual(){}		
		virtual void infoLineHeader() const;
		virtual void infoLine() const;
		virtual void refine( std::vector<LoopSet>& loopSet ); ///< Refine all clusterReps **separarly**	
		virtual void setup( LoopSet& loop ){} ///< Some derived classes may need to setup for each LoopSet
		virtual void cleanup( LoopSet& loop ){} ///< Which may then require cleanup
		void refine( LoopSet& loop );
		virtual double refine( LoopSet& loop, size_t _ClusterRepIndex ) = 0; ///< Must return the energy of the refined state
	};

	class PD_API ArcusRefine_CGMin : public ArcusRefineIndividual
	{
	public:
		ArcusRefine_CGMin() : UpdateNList(1), m_RotLib(NULL) {}

		int UpdateNList;

		inline void enableRotamerPack( const Library::RotamerLibrary& _RotLib ) { m_RotLib = &_RotLib; }
		inline void disableRotamerPack() { m_RotLib = NULL; }

	protected:
		virtual void setup( LoopSet& loop );
		virtual void cleanup( LoopSet& loop );
		virtual double refine( LoopSet& loop, size_t _ClusterRepIndex );

		const Library::RotamerLibrary* m_RotLib;

		CloneHolder< Manipulator::RotamerApplicator_SCWRL > m_SCWRL;
		CloneHolder< Physics::FF_Restraint_Positional > m_Restraint;
		CloneHolder< Protocol::Minimisation > m_Min;
	};

	class PD_API ArcusRefine_ToTra : public ArcusRefineIndividual
	{
	public:
		ArcusRefine_ToTra( const std::string& _FileStem ) : m_FileStem(_FileStem) {}

	protected:
		std::string m_FileStem;
		virtual void setup( LoopSet& loop );
		virtual double refine( LoopSet& loop, size_t _ClusterRepIndex );

		CloneHolder< IO::OutTra_BTF > m_Tra;
	};
}

#endif

