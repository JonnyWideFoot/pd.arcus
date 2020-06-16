#ifndef __FF_SOFT_STERIC_H
#define __FF_SOFT_STERIC_H

#include "forcefields/forcefield.h"

namespace Physics 
{





//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API FF_SoftVDW : public ForcefieldBase
	{
	public:
		FF_SoftVDW( WorkSpace &newwspace ):  
			ForcefieldBase( newwspace ), 
			m_Hard(false)
		{
			settodefault();
			name = "FF_SoftVDW";
		}

		virtual FF_SoftVDW* clone() const { return new FF_SoftVDW(*this); }

		virtual void settodefault()
		{
			AutoScaling14 = true;
			Hardness = 12.0;
			VdwCutoff = 15;
			VdwInnerCutoff = 12;
			Vdw14Scaling = 1.0;
			OlapOffset = 0.6;
		};

		/// How much of a VDW overlap do we allow - default 0.5 A - truncate the repulsive potential
		/// This is useful as we have no attractive VDW term, and otherwise atoms get repelled too far.
		double OlapOffset;

		/// How hard do we want our spheres?
		double Hardness;

		/// Obtain 1-4 scaling parameters from forcefield parameter set (ffps) or use defaults 
		bool AutoScaling14;

		/// Interaction cutoff distance for VdW energies [Angstrom]
		double VdwCutoff;

		/// Inner cutoff distance for VdW energies [Angstrom]
		double VdwInnerCutoff;

		/// Factor for scaling 1-4 VdW interactions (many forcefields use factor < 1.0)
		double Vdw14Scaling;

		virtual void setup();
		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void info() const; // prints a little block of parameter information
		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function

	protected:
		void calcForces_NormalVDW();
		void calcForces_Soft();

		bool m_Hard; ///< Are we a soft or hard FF?
	};
} // namespace Physics

#endif

