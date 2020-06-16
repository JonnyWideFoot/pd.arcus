#ifndef __MINIMISE_H
#define __MINIMISE_H

#include "protocols/protocolbase.h" // Provides a base class
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Protocol
{
	//-------------------------------------------------
	//
	/// \brief General Gradient based Minimisation 
	///
	/// \details Implements 2 types of Minimisation: 
	///      Cartesian Minimisation: Steepest descent and ConjugateGradient
	///
	/// \author Mike Tyka 
	///
	class PD_API Minimisation: public PickedProtocolBase 
	{
	public:
		Minimisation(Physics::Forcefield & _ff);
		Minimisation( Physics::Forcefield& _ff, const PickBase& _Picker );
		virtual Minimisation* clone() const;

		/// Run the protocol
		virtual int runcore();

		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

		/// Type of Minimisation
		enum MinType
		{ 
			SteepestDescent, 
			ConjugateGradients 
		} Algorithm;

		/// The initial stepsize. An arbitrary sort of number here, 0.1-10 conservative, 10-100 midrange, 100-10000 confident
		double  StepSize;

		/// Gradient cutoff (Energy change per step), negative means no slope cutoff
		double  SlopeCutoff;

		/// Be strict about energy going down on every step (Strictness < 0) or allow 
		/// occasional small up energy steps (Strictness > 0 )

		double  Strictness;
		
		/// Current step number
		int     Step;

	protected:
		virtual void settodefault();
	private:
		int doSteepestDescentStep();
		int doConjugateGradientStep();
		
		/// prints a little block of parameter information
		virtual void info() const;           

		/// prints a line of current energies
		virtual void infoLine() const;       


		double m_StepMultiplier;
		double m_OldEnergy;
	};

} // namespace 'Protocol'

#endif

