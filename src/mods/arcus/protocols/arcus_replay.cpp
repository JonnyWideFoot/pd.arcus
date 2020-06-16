#include "global.h"
#include "arcus_replay.h"

namespace Protocol
{
	ArcusReplay::ArcusReplay( Physics::Forcefield& _ffs, Physics::Forcefield& _ff, const Library::AngleSet& _as, const Library::RotamerLibrary& _rotLib, InputTrajectory& tra )
		: ArcusBase( _ffs, _ff, _as ), 
		m_RotLib(&_rotLib),
		m_Tra(&tra)
	{
		setDefaultRefiner();
		m_SplitBuilder = true; // most likely the best plan
	}

	void ArcusReplay::setDefaultRefiner()
	{
		m_Stage3Refiner = &m_DefaultStage3;
	}

	void ArcusReplay::setRefiner( ArcusRefineBase& _refiner )
	{
		m_Stage3Refiner = &_refiner;
	}

	int ArcusReplay::runcore()
	{
		WorkSpace& wspace = getWSpace();
		SnapShot snap;
		m_Tra->reset();

		// ------------------------------------------
		//  Part 1: Create the m_BuiltSections array
		// ------------------------------------------

		std::vector<PosStore> posCache;
		for( size_t q = 0; q < m_Regions.size(); q++ )
		{
			PosStore curentPos( getWSpace(), m_Regions[q] );
			posCache.push_back( curentPos );
		}

		m_BuiltSections.clear();
		m_BuiltSections.resize(m_Regions.size());
		for( size_t q = 0; q < m_Regions.size(); q++ )
		{
			LoopSet& loopSet = m_BuiltSections[q];
			loopSet.range = PickResidueRange( m_Regions[q] );
		}

		int modelCount = 0;
		while( m_Tra->readNext( snap ) )
		{
			modelCount++;
			wspace.load( snap );
			for( size_t q = 0; q < m_Regions.size(); q++ )
			{			
				posCache[q].store();
				LoopSet& loopSet = m_BuiltSections[q];	
				loopSet.posCache.push_back(posCache[q]);				
			}
		}

		// Set pointers **AFTER** vector allocation above! Otherwise pointers could be invalid!!
		for( size_t q = 0; q < m_Regions.size(); q++ )
		{
			LoopSet& loopSet = m_BuiltSections[q];	
			for( size_t r = 0; r < loopSet.posCache.size(); r++ )
			{		
				loopSet.clusterReps.push_back( r ); // just a sequential index here
				loopSet.posCachePointer.push_back( 
					std::pair<double,PosStore*>( 0.0, &loopSet.posCache[ r ] ));
			}
		}

		// ---------------------------------------------------------------
		//  Part 2: Now call the refiner on the m_BuiltSections array :-D
		// ---------------------------------------------------------------

		wspace.Step = 0; // this needs to be reset
		m_Stage3Refiner->enableComments( getTraComments() );
		m_Stage3Refiner->initialise( getWSpace(), *ff );
		m_Stage3Refiner->refine( m_BuiltSections ); // Stage-3 (bridge-pattern)

		return modelCount;
	}
}

