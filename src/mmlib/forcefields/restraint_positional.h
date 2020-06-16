#ifndef __FFRESTRAINT_POSITIONAL_H
#define __FFRESTRAINT_POSITIONAL_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class

namespace Physics
{

	//-------------------------------------------------
	//
	/// \brief Implements a simple Cartesian harmonic, positional restraint for sets of atoms
	///
	/// \details 
	///    
	///
	/// \author Mike Tyka  
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API FF_Restraint_Positional : public AtomRestraintForcefieldBase
	{
	public:
		FF_Restraint_Positional(WorkSpace &newwspace)
			: AtomRestraintForcefieldBase(newwspace) 
		  {
				name = "Positional Rest.";
				ShortName = "PosRest";

			  k = 0;
			  Power = 2;
		  };

		virtual FF_Restraint_Positional* clone() const { return new FF_Restraint_Positional(*this); }


		/// Power of the restraint (epot = k*|x-x0|^Power )
		double Power; 

		virtual double calcEnergyAtQ(double newQ);

		/// prints a little block of parameter information
		virtual void info() const; 

	protected:
		virtual void setup();

		virtual void  calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void  calcEnergies();
		virtual void  calcForces();

	};



} // namespace Physics


#endif

