#ifndef __FFRESTRAINT_INTERNAL_H
#define __FFRESTRAINT_INTERNAL_H

#include "forcefields/forcefield.h" // provides a base class
#include "forcefields/restraintbase.h" // provides a base class
#include "forcefields/ffbonded.h"   // provides a member variable

namespace Physics
{


	//-------------------------------------------------
	//
	/// \brief Restraint to maintain a structure by the addition of internal harmonic struts within a set of atoms 
	///
	/// \details This Forcefield adds internal struts which keep a structure rigid internally, without affecting
	///          the rigid body movement of the set of atoms selected. Struts are added to all pairs of atoms 
	///          within the selected set of atoms which are within the parameter RestCutoff.  
	///
	/// \author Mike Tyka  
	///
	/// \todo 
	///
	/// \bug 
	///
	class PD_API FF_Restraint_Internal : public AtomRestraintForcefieldBase
	{
	public:
		FF_Restraint_Internal( WorkSpace &newwspace );
		virtual FF_Restraint_Internal* clone() const { return new FF_Restraint_Internal(*this); }

		virtual void info() const; // prints a little block of parameter information

		double RestCutoff;

		bool DivByNumber;       // divide k by number of pairs ??

	private:
		virtual void setup();

		int createRestraints();

		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();
		virtual double calcEnergyAtQ(double newQ);

		std::vector< Bond > rest;
		int nrest; // number of restraints
	};


} // namespace Physics


#endif

