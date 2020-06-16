#ifndef __FFRESTRAINT_TORSIONAL_H
#define __FFRESTRAINT_TORSIONAL_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class

#include "forcefields/ffbonded.h"   // provides a member variable

namespace Physics
{
	class PD_API Torsion;

	//-------------------------------------------------
	//
	/// \brief Restrains all or some of the internal torsions within a molecule 
	///
	/// \details 
	///    
	/// \author Mike Tyka  
	///
	/// \todo HALFBAKED 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_Torsional : public AtomRestraintForcefieldBase
	{
	public:
		FF_Restraint_Torsional(  WorkSpace &newwspace );
		virtual FF_Restraint_Torsional* clone() const { return new FF_Restraint_Torsional(*this); }

		virtual double calcEnergyAtQ(double newQ);
		
		/// prints a little block of parameter information
		virtual void info(); 
		
		/// prints a little block of parameter information
		virtual void detail(); 

		bool OneRestraintPerBond; 
	protected:
		virtual void setup();
		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

	private:
	 std::vector< Torsion > torsion;
	};


} // namespace Physics


#endif

