#ifndef __ARCUS_REPLAY_H
#define __ARCUS_REPLAY_H

#include "arcus_refiner.h"

namespace Library
{
	class RotamerLibrary;
}

class InputTrajectory;

namespace Protocol
{
	class PD_API ArcusReplay : public ArcusBase, public IO::BTF_Block_Comment_User
	{
	public:
		ArcusReplay( 
			Physics::Forcefield& _ffs, 
			Physics::Forcefield& _ff, 
			const Library::AngleSet& _as, const 
			Library::RotamerLibrary& _rotLib,
			InputTrajectory& _tra );

		virtual ArcusReplay* clone() const 
		{ 
			return new ArcusReplay(*this); 
		}

		virtual int runcore();

		/// Do nothing
		virtual int initialise()
		{
			return ArcusBase::initialise();
		}

		virtual void configLevel1Filters( size_t i )
		{
			// Allocate it, but add nothing, we have no stage-1
			m_Filters.push_back( FilterContainer() );
		}

		void setDefaultRefiner();
		void setRefiner( ArcusRefineBase& _refiner );

	protected:
		const Library::RotamerLibrary* m_RotLib;
		InputTrajectory* m_Tra;
		ArcusRefineBase* m_Stage3Refiner;
		ArcusRefine_CGMin m_DefaultStage3;
		std::vector< LoopSet > m_BuiltSections;
	};
}

#endif

