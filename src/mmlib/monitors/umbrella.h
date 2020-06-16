#ifndef __UMBRELLA_H
#define __UMBRELLA_H

#include "monitors/monitorbase.h" // Provides base class
#include "forcefields/forcefield.h"
#include "forcefields/restraintbase.h"

namespace Monitors{
	class PD_API UmbrellaMonitor;

	void printUmbrellaProfile(UmbrellaMonitor &umbrella, double binwidth);






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API UmbrellaMonitor: public MonitorBase{
	public:

		// SWIG seems to get confused here ..
#ifndef SWIG
		friend void printUmbrellaProfile(UmbrellaMonitor &umbrella, double binwidth);
#endif

		UmbrellaMonitor(Physics::RestraintForcefieldBase *newff )
		{
			name = "umbr";
			ff = newff;
		}

		void setcurdata(){
			addData(  ff->Q );
		}
	private:
		Physics::RestraintForcefieldBase *ff;
	};
}
#endif


