#include "global.h"
#include "library/angleset.h"
#include "workspace/workspace.h" 
#include "workspace/neighbourlist.h"
#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "protocols/torsionalminimisation.h"
#include "protocols/minimise.h"
#include "arcus_stitch.h"

namespace Protocol
{
	ArcusStitch::ArcusStitch( Physics::Forcefield& _ffs, Physics::Forcefield& _ff, Library::AngleSet& _as )
		: ArcusBase( _ffs, _ff, _as )
	{
	}

	int ArcusStitch::runcore()
	{
		// Print a header line or some other information
		info();
		infoLineHeader();

		// What are we actually moving?
		PickAtomRanges picker;
		setPickedRegions( getWSpace(), picker, m_Regions );

		// Energy assessment / best structure keep
		double bestEpot = DBL_MAX;
		PosStore bestFound(getWSpace(),picker);
		bestFound.store();

		// Main simulation loop, one model per step
		for(Step = 0; Step < Steps; Step++) 
		{
			if( OutputLevel >= Verbosity::Normal )
				Printf("Conformer itteration %d enumeration:\n")(Step);

			// Set initial valid conformers
			for( size_t j = 0; j < m_Regions.size(); j++ )
			{
				// Manipulate each section in turn
				size_t conformerAttempts = 0;
				while( m_ConfBuilder[j].next() )
				{
					conformerAttempts++;
					this->getWSpace().outtra.append();
					if( m_Filters[j].passes() )
					{
						break;
					}
				}
				m_ConfBuilder[j].reset();
				if( OutputLevel >= Verbosity::Normal )
					Printf("  Section %d required %d attempts to pass filters\n")(j)(conformerAttempts);
			}

			// Fast torsional under steric FF on all regions sumultaneously!
			// Use default parameters...
			Protocol::TorsionalMinimisation torMin( *ffs, picker );
			torMin.Steps = 500;
			torMin.SlopeCutoff = 0.1 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
			torMin.InitialCapFactor = 0.1;
			torMin.UpdateScr = 10;
			torMin.UpdateTra = 10;
			torMin.run();

			// Slower cartesian on big ff
			Minimisation min2( *ff, picker );
			min2.Steps = 1501; // *maximum*, but will be less than this as 'SlopeCutoff' is defined
			min2.SlopeCutoff = 0.01 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
			min2.StepSize = 2E1;
			min2.UpdateScr = 10;
			min2.UpdateTra = 10;
			min2.run();
			if( min2.Step == min2.Steps ) 
				printf("WARNING: Minimisation maxed out the Step allocation\n");

			if( getWSpace().ene.epot < bestEpot )
			{
				bestFound.store();
				bestEpot = getWSpace().ene.epot;
			}

			// Call any monitors that might have been added to you as a protocol
			runmonitors();

			// Display your infoline every so often (UpdateScr)
			if((OutputLevel)&&(UpdateScr > 0) && (((Step) % UpdateScr) == 0))
			{
				infoLine();
			}			

			// Save the coordinates in the trajectory as required
			if((OutputLevel)&&(UpdateTra > 0)) 
			{
				getWSpace().outtra.append();
			}		
		}

		// print some statistics if you want
		if(OutputLevel) 
		{
			printFinalStatistics(Steps);
		}

		if( FinalState == LowestEpot )
		{
			bestFound.revert();
		}
		// else its LastAcc, which we have already...

		// runcore should return the number of force/energy evaluations 
		return Step;
	}
}

