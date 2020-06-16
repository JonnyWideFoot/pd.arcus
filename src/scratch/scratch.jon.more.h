#ifndef __SCRACTH_JON_MORE_H
#define __SCRACTH_JON_MORE_H

#include "global.h"

#include "forcefields/forcefield.h"
#include "protocols/montecarlo.h"
#include "workspace/workspace.fwd.h"

#include "scratch.jon.basics.h"

namespace Protocol
{
	class MonteCarloLowKeep : public MonteCarlo
	{
	public:
		MonteCarloLowKeep(ProtocolBase&  _evaluator,
		Manipulator::MoveSet& _moveset);
		virtual int runcore();

	protected:
		// Redefine acceptance function
		virtual bool accept(
			double enenew,
			double eneold
		) const;

	private:
		double m_StartEne;
		mutable double m_BestEne;
	};
}

void Test_Conformer();
//void TestSoftSteric( WorkSpace& wspace );
//void TestSoftSteric();
//void Test_CoiledCoil();
//void TestPops();
//void Test_Mike();
//void main_CraigSim();
//int Tweaky(int argc, char** argv);
void main_GentleHarmonicRestraintMinimisation();
//void TempFunc();
void TestProximityGrid();
void randomRotamerApplication( bool _testSterics, bool importForeignLib );
void TheWandersOfRotamers();
void testGraphTheory();
void doMD();
void coiledCoilMake();
void zincFingerPrimaryModelBuilder();

#endif

