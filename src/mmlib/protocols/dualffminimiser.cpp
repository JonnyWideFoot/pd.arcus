#include "global.h"
#include "forcefields/forcefield.h"
#include "minimise.h"
#include "dualffminimiser.h"

using namespace Physics;

namespace Protocol
{
	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff ) 
		: Minimisation( _ff ), stericff(NULL)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff)  
		: Minimisation( _ff ), stericff(&_stericff)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, const PickBase& _Picker )  
		: Minimisation( _ff, _Picker ), stericff(NULL)
	{
		setDefaults();
	}

	DualFFMinimiser::DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff, const PickBase& _Picker )
		: Minimisation( _ff, _Picker ), stericff(&_stericff)
	{
		setDefaults();
	}

	void DualFFMinimiser::setDefaults()
	{
		SDPreMinSteps = 0;
		StericMinSteps = 0;
		StericStepSize = 2E5;
		StericSlopeCutoff = -1.0;
		StericKillFull = DBL_MAX;
	}

	int DualFFMinimiser::runcore()
	{
		double unitConvStericKillFull = StericKillFull * (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na);

		int Steps = 0;

		if( SDPreMinSteps > 0 )
		{
			// Perform a steepest descent preminimisation using the steric ff if 
			// available and the normal ff if not.
			Physics::Forcefield* ffSend = (stericff!=NULL) ? stericff : ff;
			ASSERT( ffSend != NULL, CodeException, "DualFFMinimiser internal NULL reference!");
			if( OutputLevel ) printf("Steep Pre-Minimisation:\n");
			Minimisation prepreminimise( *ffSend, getPicker() );
			prepreminimise.Algorithm = Minimisation::SteepestDescent;
			prepreminimise.Steps = SDPreMinSteps;

			// Mirror parent update params
			prepreminimise.UpdateScr = UpdateScr;
			prepreminimise.UpdateTra = UpdateTra;
			prepreminimise.UpdateMon = UpdateMon;
			prepreminimise.UpdateNList = UpdateNList;
			prepreminimise.OutputLevel = OutputLevel;

			Steps += prepreminimise.run();
		}

		if( StericMinSteps > 0 )
		{
			ASSERT( stericff != NULL, ArgumentException, "StericMinSteps can only be used if a steric forcefield has been provided");
			if( OutputLevel ) printf("Steric Minimisation:\n");
			Minimisation stericMinimise( *stericff, getPicker() );
			stericMinimise.Algorithm = Minimisation::ConjugateGradients;
			stericMinimise.Steps = StericMinSteps;
			stericMinimise.StepSize = StericStepSize;
			stericMinimise.SlopeCutoff = StericSlopeCutoff;

			// Mirror parent update params
			stericMinimise.UpdateScr = UpdateScr;
			stericMinimise.UpdateTra = UpdateTra;
			stericMinimise.UpdateMon = UpdateMon;
			stericMinimise.UpdateNList = UpdateNList;
			stericMinimise.OutputLevel = OutputLevel;

			Steps += stericMinimise.run();
		}
		else if( stericff != NULL )
		{
			printf("CODE WARNING: DualFFMinimiser stericff is defined, but StericMinSteps is '0'. No steric minimisation has been performed!\n");
		}

		ff->calcEnergies(); // trigger a single energy calc on the full forcefield
		// Then test to see we have done what we have resolved steric issues
		if( getWSpace().ene.epot < unitConvStericKillFull )
		{
			if( OutputLevel ) 
			{
				printf("Full Minimisation:\n");
			}
			Steps += Minimisation::runcore();
		}
		else
		{
			if( OutputLevel ) 				
			{
				printf("Full Minimisation: Killed as Steric did not succeed. Ene:'%8.3f', Cutoff:'%8.3f'\n",
					getWSpace().ene.epot * Physics::PhysicsConst::Na / Physics::PhysicsConst::kcal2J, StericKillFull);
			}
		}

		return Steps;
	}
} // namespace




