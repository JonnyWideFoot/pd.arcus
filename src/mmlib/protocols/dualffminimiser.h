// -------------------------------------------------------------------------------------
// Description: A minimisation protocol that is very tollerant of initial steric clash.
// --------------------------------------------------------------------------------------

#ifndef __DUAL_FF_MINIMISER_H
#define __DUAL_FF_MINIMISER_H

// Essential Headers
#include "protocols/minimise.h"
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"


namespace Protocol
{





//-------------------------------------------------
//
/// \brief A combine protocol which does two minimisations on different forcefields 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
///

	class PD_API DualFFMinimiser: public Minimisation
	{
	public:
		
		DualFFMinimiser( Physics::Forcefield & _ff );
		DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff);
		DualFFMinimiser( Physics::Forcefield & _ff, const PickBase& _Picker );
		DualFFMinimiser( Physics::Forcefield & _ff, Physics::Forcefield & _stericff, const PickBase& _Picker );

		virtual DualFFMinimiser* clone() const 
		{ 
			return new DualFFMinimiser(*this); 
		}


		virtual int runcore();


		/// Perform a quick Steepest Descent Minimisation first?
		int SDPreMinSteps; 

		/// Perform a steric conjugate gradients Minimisation?
		int StericMinSteps; 

		/// If the steric minimisation is enabled and doesnt get below this energy, the full minim is not performed.
		double StericKillFull; 

		/// The size of the step size in the steric forcefield minimisation
		double StericStepSize; 

		/// The energy gradient cutoff uses in the steric forcefield minimisation
		double StericSlopeCutoff; 

	protected:
		void setDefaults();
	private:
		Physics::Forcefield* stericff; // steric forcefield if one was given
	};
}

#endif

