// --------------------------------------------------------------------------------
// File: csa.h
// Description: Conformational Space Annealing (Scheraga et al. 1999)
// Version: 0.1
// Author: Michael Tyka
// mail: m.tyka@bristol.ac.uk
//
// Dependencies: wspace.h
// Headers: csa.h
// Sources: csa.cpp
// Classe(s): ConformationalSpaceAnnealing
//

#ifndef __CSA_H
#define __CSA_H

// Essential Headers
#include "protocols/protocolbase.h"
#include "forcefields/forcefield.h"
#include "manipulators/basicmoves.h"

// Forward Declarations
namespace Library
{
	class PD_API RotamerLibrary;
	class PD_API AngleSet;
}

class PD_API SnapShot;


namespace Protocol{

	class ConformationalSpaceAnnealingPS
	{
	public:
		ConformationalSpaceAnnealingPS(){
			settodefault();
		};

		int UpdateScr;
		int UpdateTra;
		int Steps;
		double Temperature;
		int scramble;
		int librarysize;
		int maxgroupsize;
		int repeats;
		int mcmsteps;
		double mcmTemperature;
		int mcmpreminsteps;
		int mcmminsteps;
		double rmslimit;
		double rmslimit_end;
		double exchange_prop;
		int finalsasteps; // do a simulated annealing on the final librayr?
		double finalsastartt; // start and end Temperatures ?
		double finalsaendt;
		int finalsaminsteps; // and a minimisation ?
		bool optimiseonrepeat1;
		std::string stem; // filename stem for csa files

		virtual void settodefault(){
			Steps = 0;
			UpdateScr = 0;
			UpdateTra = 0;
			Temperature = 300;
			scramble = 1;
			librarysize = 0;
			maxgroupsize = 10;
			repeats = 1;
			mcmsteps = 1;
			mcmTemperature = 300;
			mcmpreminsteps = 0;
			mcmminsteps = 1;
			rmslimit = 2.0;
			rmslimit_end = 2.0;
			exchange_prop = 0.2;
			finalsasteps = 0; // do a simulated annealing on the final librayr?
			finalsastartt = 300; // start and end Temperatures ?
			finalsaendt = 100;
			finalsaminsteps = 1000; // and a minimisation ?
			stem = "csastem";
			optimiseonrepeat1 = false;
		};
		void info();
	private:
	};

	class ConformationalSpaceAnnealing: public ProtocolBase
	{
	public:
		ConformationalSpaceAnnealing(WorkSpace * wspace,
			Physics::Forcefield * _ff,
			Physics::Forcefield * _stericff,
			Library::AngleSet * _angset,
			std::vector < Manipulator::MoveBase * >* _sp_csa,
			std::vector < Manipulator::MoveBase * >* _sp_fine,
			ConformationalSpaceAnnealingPS &_ps);

		~ConformationalSpaceAnnealing();

		virtual ConformationalSpaceAnnealing *clone()const{
			return new ConformationalSpaceAnnealing(*this);
		}

		virtual int runcore();

		virtual void info() const; // prints a little block of parameter information
		virtual void infoLine() const; // prints a line of current energies
		virtual void infoLineHeader() const; // prints the headers for the above function

		void setFinalSAParams(
			int _finalsasteps, // do a simulated annealing on the final librayr?
			double _finalsastartt, // start and end Temperatures ?
			double _finalsaendt,
			int _finalsaminsteps ){ // and a minimisation ?
				ps.finalsasteps = _finalsasteps;
				ps.finalsastartt = _finalsastartt;
				ps.finalsaendt = _finalsaendt;
				ps.finalsaminsteps = _finalsaminsteps;
		}

		int fillLibrary(Library::RotamerLibrary *rotlib); // fills library with current wspace structure;

	private:

		Protocol::ProtocolBase * mimimiser;
		Physics::Forcefield * stericff;
		std::vector < Manipulator::MoveBase * >* sp_csa;
		std::vector < Manipulator::MoveBase * >* sp_fine;

		Library::AngleSet* angset;

		int Step;

		SnapShot *psplib;

		ConformationalSpaceAnnealingPS ps; // CSA parameter set
	};

} // namespace 'Protocol'

#endif
