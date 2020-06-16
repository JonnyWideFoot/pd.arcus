#include "global.h"

#include "manipulators/movebase.h"
#include "protocols/temperature.h"
#include "workspace/workspace.h"
#include "forcefields/forcefield.h"

#include "montecarlo.h"

using namespace Maths;
using namespace Physics;
using namespace Manipulator;

namespace Protocol
{
	int MonteCarlo::runcore()
	{
		if(OutputLevel)
		{
			if((UpdateScr>0)&&(UpdateScr<Steps))
			{
				infoLineHeader();
				ff->infoLineHeader();
			}
		}

		bool doReasons = OutputLevel && ReportFilterFailReasons;

		reset();

		// save original state
		oldstate = getWSpace().save();
		lowstate = oldstate;	

		for(Step = 0; Step < Steps; Step++) 
		{
			getWSpace().Step = Step;
			if(every(Step, 100)) getWSpace().cleanSpace(); // occasionally move any stray molecules back into simulation

			// make a change to the structure
			getWSpace().resetAllAtomMovedFlags();
			if(Step != 0) moveset->apply();

			if( DoZeroGeometry ) getWSpace().zeroCentreOfGeometry();

			validity = Reject; // RESET to INVALID

			if( PreFilters.passes() ) // PRE-screen PRIOR to using the evaluator (thats the main point of a filter)
			{		
				evaluations += evaluator->runcore();
				// reject/accept structure according to metropolis criterion
				bool accepted = accept(getWSpace().ene.epot, oldstate.epot);

				// If we have already failed, there is no point in running post-filters.
				if( accepted )
				{
					// Do we also pass any post-filters?
					accepted = PostFilters.passes();
					if( !accepted )
					{
						// Well wasn't that a waste of an evaluation, maybe we should pre-screen better ;-)
						validity = PostFilterFail; 
						if( doReasons ) 
						{ 
							std::cout << PostFilters.reason() << std::endl; 
						}
					}										
					else if( lowestEne == DBL_MAX || lowestEne > getWSpace().ene.epot )
					{
						lowstate = getWSpace().save();
						lowestEne = getWSpace().ene.epot;
						validity = Accept; // Cool!
					}
					else
					{
						validity = Accept; // Decent
					}
				}
			}
			else
			{
				validity = PreFilterFail;
				if( doReasons ) 
				{
					std::cout << PreFilters.reason() << std::endl;
				}
			}

			// Cosmetics:
			if((OutputLevel) && every(Step,UpdateScr)) 
			{
				infoLine();
			}

			if(validity != Accept)
			{				
				getWSpace().load(oldstate);
				nonacceptances++;
				repeatmonitors();
				if(UpdateTraRej && every(Step,UpdateTra)) getWSpace().outtra.append();
			}
			else
			{
				
				oldstate = getWSpace().save();
				acceptances++;
				nonacceptances = 0;
				runmonitors(); // update monitors with a new measurement 
				if(UpdateTraAcc && every(Step,UpdateTra)) getWSpace().outtra.append();
			}
		}

		switch(FinalState)
		{
		case MonteCarlo::LastAcc:
			{
				getWSpace().load(oldstate);
				break;
			}
		case MonteCarlo::LowestEpot:
			{
				if(OutputLevel) 
				{
					printf("Lowest ene: %lf \n", lowstate.epot * PhysicsConst::J2kcal * PhysicsConst::Na);
				}
				getWSpace().load(lowstate);
				break;
			}
		case MonteCarlo::Last:
			{
				break;
			}
		}

		if(OutputLevel){
			printf("mc.stats:  Steps           : %7d \n",Steps);
			printf("mc.stats:  Acceptance Rate : %8.5f \n", (double) acceptances / double (Step + 1) );
		}

		return Step;
	}

	bool MonteCarlo::accept( double enenew, double eneold ) const
	{ 
		if(!isNumber(enenew)) return false;
		if(enenew < eneold)
		{ 
			return true;
		}
		else 
		{
			double enediff = (eneold - enenew) / 
				(PhysicsConst::kB * 
				Temperature->get(double(Step)/double(Steps)));
			if(enediff < -50.0) return false;
			double expEdivkT = exp(enediff);
			return (frand() < expEdivkT);
		}
	}

	void MonteCarlo::info() const
	{
		printf("montecarlo.Steps     %d",Steps);
	}

	void MonteCarlo::addPreFilter(FilterBase &filter)
	{
		PreFilters.add(filter);
	}

	void MonteCarlo::addPostFilter(FilterBase &filter)
	{
		PostFilters.add(filter);
	}

	void MonteCarlo::infoLineHeader() const 
	{ 
		printf("Type Step evals  accs accScope epot oldepot lowepot accfreq temp\n");
	}

	void MonteCarlo::infoLine() const 
	{
		if(!OutputLevel) return;

		char accType;
		switch( validity )
		{
		case Accept:
			if( UpdateScr != 1 && !UpdateScrAcc ) return;
			accType = 'A';			
			break;
		case Reject:
			if( UpdateScr != 1 && !UpdateScrRej ) return;
			accType = 'R';
			break;
		case PreFilterFail:
			if( UpdateScr != 1 && !UpdateScrRej ) return;
			accType = '<';
			break;
		case PostFilterFail:
			if( UpdateScr != 1 && !UpdateScrRej ) return;
			accType = '>';
			break;
		default:
			THROW(CodeException,"Unknown internal validicy class in MonteCarlo");
		}

		double accScope = (double) (acceptances-acceptances_block) / double (UpdateScr);

		printf("%c%6d%8d%4d %6.2lf %9.3lf %9.3lf %9.3lf %5.1lf %5.1lf",
			accType, Step, evaluations, acceptances, accScope,
			getWSpace().ene.epot * PhysicsConst::J2kcal * PhysicsConst::Na,
			oldstate.epot * PhysicsConst::J2kcal * PhysicsConst::Na,
			lowestEne * PhysicsConst::J2kcal * PhysicsConst::Na,
			(double) acceptances / double (Step + 1),
			Temperature->get(double(Step)/double(Steps)));
		mon.printCurData();
		ff->infoLine();
		printf("\n");

		// Update the block over the last UpdateScr steps
		acceptances_block = acceptances;
	}
} // namespace 'Protocol'

