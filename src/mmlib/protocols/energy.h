#ifndef __ENERGY_H
#define __ENERGY_H

// Essential Headers
#include "protocols/protocolbase.h" // Provides a base class
#include "workspace/snapshot.h" // Provides a class member
#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Protocol
{

//-------------------------------------------------
//
/// \brief  A passive protocol that calculates the energy 
///
/// \details  runcore does nothing except for calling calcEnergy of the supplied Forcefield
///
/// \author Mike Tyka 
///
	class PD_API Energy: public ProtocolBase {
	public:
		Energy(	Physics::Forcefield & _ff):
		ProtocolBase(_ff)
		{ };

		virtual Energy* clone() const 
		{ 
			return new Energy(*this); 
		}

		virtual ~Energy(){};

		/// runs quietly, calculates energy and leaves
		virtual int runcore(); 

	private:

		/// prints a little block of parameter information
		virtual void info() const{}; 

		/// prints a line of current energies
		virtual void infoLine() const{}; 

		/// prints the headers for the above function
		virtual void infoLineHeader() const{}; 
	};


} // namespace 'Protocol'

#endif

