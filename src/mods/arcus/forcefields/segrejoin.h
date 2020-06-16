#ifndef __LOOP_JOIN_ASSESSMENT_H
#define __LOOP_JOIN_ASSESSMENT_H

// mmlib
#include "forcefields/restraintbase.h"
#include "forcefields/ffcustom.h"
#include "filters/basicfilters.h"

// loopbuilder
#include "workspace/segdef.h"

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
	class PeptideGroupRejoinForce : public SegmentDefUser, public ForcefieldBase
	{
	public:
		PeptideGroupRejoinForce(SegmentDef& _segdef, Verbosity::Type _verbosity = Verbosity::Normal );
		virtual PeptideGroupRejoinForce* clone() const { return new PeptideGroupRejoinForce(*this); }
		bool Initialise( double strength = 50.0, bool useCARestraints = false );

		virtual void info() const { m_ForceField.info(); } // prints a little block of parameter information
		virtual void infoLine() const { m_ForceField.infoLine(); } // prints a line of current energies
		virtual void infoLineHeader() const { m_ForceField.infoLineHeader(); } // prints the headers for the above function

		virtual void calcForces();

		void enableHarmonicCARestraints( bool enabled );

	private:
		void addHarmonicCARestraints();
		void updateHarmonicCARestraints();

		virtual void setup(); // obtain the internal atom pointers and make the internal forcefield
		Physics::FF_Custom m_ForceField;
		bool m_CARestraintsOn;
		double m_RestraintStrength; // seems empiricly good. Too strong and it wont search, too weak and it wont get there.
	};
}

#endif


