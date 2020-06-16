#include "global.h"
#include "workspace/workspace.h"
#include "forcefields/forcefield.h"
#include "library/rotamerlib.h"
#include "manipulators/movebase.h"
#include "manipulators/basicmoves.h"
#include "protocols/dualffminimiser.h"
#include "protocols/montecarlo.h"
#include "scpack.h"

void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const Library::RotamerLibrary& rotLib )
{
	// Our moveset
	Manipulator::MoveSet moves( wspace );
	Manipulator::SidechainRotamerLibMove* rotamer = new Manipulator::SidechainRotamerLibMove( wspace, rotLib, 1.0, 1.0 );
	moves.addWithOwnership( rotamer );

	// Minimise only the sidechains
	Protocol::DualFFMinimiser eval(ff, ffs, PickSidechains() );
	eval.SDPreMinSteps = 31;
	eval.StericMinSteps = 501;
	eval.Steps = 501;
	eval.StepSize = 2E1;
	//eval.!OutputLevel = true;
	eval.UpdateScr = 50;
	//eval.UpdateTra = 1;
	//eval.UpdateMon = 10;
	eval.run();

	// Perform a montecarlo procedure
	Protocol::MonteCarlo mc( eval, moves );
	mc.Steps = 10;
	mc.UpdateTra = 1;
	mc.UpdateScr = 1;
	mc.UpdateTraAcc = true;
	//mc.UpdateTraRej = true;
	mc.FinalState = Protocol::MonteCarlo::LowestEpot;
	mc.run(); // Pack sidechains using the rotamers	
}

void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const std::string& rotLibPath )
{
	Library::RotamerLibrary rotLib(wspace.ffps());
	rotLib.readLib( rotLibPath );
	MCPackSideChains( wspace, ffs, ff, rotLib );
}

