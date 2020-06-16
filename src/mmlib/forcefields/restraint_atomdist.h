#ifndef __FFRESTRAINT_ATOMDIST_H
#define __FFRESTRAINT_ATOMDIST_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class
#include "forcefields/ffcustom.h"   // provides a member variable

namespace Physics
{

	//-------------------------------------------------
	//
	/// \brief Simple Harmonic Restrant between two atoms 
	///
	/// \details 
	///
	/// \author Mike Tyka  
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_AtomDistance : public RestraintForcefieldBase
	{
	public:
		FF_Restraint_AtomDistance(WorkSpace &newwspace);

		virtual FF_Restraint_AtomDistance* clone() const 
		{ 
			return new FF_Restraint_AtomDistance(*this); 
		}

		/// prints a little block of parameter information
		virtual void info(); 
		
		/// prints a little block of parameter information
		virtual void detail() { info(); }; 

		/// restraint Power (default = 2, i.e. harmonic)
		int Power;  

		/// atom index 1
		int Atom_i;  

		/// atom indices;
		int Atom_j; 

		/// equilibrium distance
		double Dist_ij; 

	private:
		virtual void setup();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();
		virtual double calcEnergyAtQ(double newQ);

		virtual void infoLine() const;
		virtual void infoLineHeader() const; // prints the headers for the above function

		FF_Custom cff;
	};


} // namespace Physics


#endif

