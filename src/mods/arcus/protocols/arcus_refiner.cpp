#include "global.h"
#include "fileio/tra.h"
#include "arcus_refiner.h"

namespace Protocol
{
	void ArcusRefineIndividual::refine( std::vector<LoopSet>& loopSet )
	{
		ASSERT( m_WSpace != NULL, NullInternalException, "ArcusRefineBase is uninitialised, WorkSpace pointer is NULL");
		ASSERT( m_FF != NULL, NullInternalException, "ArcusRefineBase is uninitialised, Forcefield pointer is NULL");

		// At the moment, we perform independent refinement ....
		for( size_t j = 0; j < loopSet.size(); j++ )
		{
			setup(loopSet[j]);
			refine(loopSet[j]);
			cleanup(loopSet[j]);
		}		
	}

	void ArcusRefineBase::initialise( WorkSpace& _wspace, Physics::Forcefield& _ff )
	{
		m_WSpace = &_wspace;
		m_FF = &_ff;
	}

	void ArcusRefineIndividual::infoLineHeader() const 
	{ 
		// prints the headers for the above function
		printf("%8s\t%6s\t%10s\t", "Step", "CRMS", "Epot");
		m_FF->infoLineHeader();
		printf("\n");
	}

	void ArcusRefineIndividual::infoLine() const 
	{ 
		// prints a line of current energies
		printf("%8d\t%6.2lf\t%15.10lf\t",
			getWSpace().Step, getWSpace().ene.cRMS, (double)getWSpace().ene.epot * 
				Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na);
		m_FF->infoLine();
		printf("\n");
	}

	void ArcusRefineIndividual::refine( LoopSet& loopSet )
	{
		double bestEne = DBL_MAX;
		PosStore cacheBest( getWSpace(), loopSet.range );

		ASSERT( loopSet.clusterReps.size() > 0, ProcedureException, "No data for ArcusRefineIndividual::refine()" );

		if( OutputLevel ) infoLineHeader();
		for( size_t j = 0; j < loopSet.clusterReps.size(); j++ )
		{
			double ene = refine( loopSet, j );
			if( ene < bestEne )
			{
				bestEne = ene;
				cacheBest.store();
			}
			if( OutputLevel ) infoLine();
		}

		cacheBest.revert();

		if( AppendBestToTra )
		{
			if( hasTraComments() ) 
				getTraComments()->setFormat( "BEST found: Minimised-ene is: %8.3lf. CRMS: %8.3lf" )
				(bestEne)
				(cacheBest.calcCRMS(PickHeavyAtoms()));				
			getWSpace().outtra.append();
		}
	}

	void ArcusRefine_CGMin::setup( LoopSet& loopSet )
	{
		// Config the Minimisation protocol
		m_Min = new Protocol::Minimisation( *m_FF, loopSet.range );
		m_Min->Steps = 1500;
		m_Min->UpdateNList = UpdateNList; // propogate the parameter to the child-protocol
		m_Min->SlopeCutoff = 0.01 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
		m_Min->OutputLevel = Verbosity::Silent;
		//m_Min->UpdateScr = 50;
		//m_Min->UpdateTra = 10;

		// Restrain cartesian
		m_Restraint = new Physics::FF_Restraint_Positional( getWSpace() );
		m_Restraint->k = 20.0;
		m_FF->add( m_Restraint.data() ); // popped during cleanup below...

		if( m_RotLib != NULL )
		{
			m_SCWRL = new Manipulator::RotamerApplicator_SCWRL( 
				getWSpace(), 
				*m_RotLib,
				PickResidueList(loopSet.range) );	
			m_SCWRL->OutputLevel = Verbosity::Silent;
		}
		else
		{
			m_SCWRL.release();
		}
	}
	
	void ArcusRefine_CGMin::cleanup( LoopSet& loopSet )
	{
		m_FF->pop_back();
	}

	double ArcusRefine_CGMin::refine( LoopSet& loopSet, size_t _ClusterRepIndex )
	{
		int outerStep = getWSpace().Step;

		size_t i = loopSet.clusterReps[_ClusterRepIndex];
		loopSet.posCachePointer[i].second->revert();

		if( m_SCWRL.assigned() )
		{
			m_SCWRL->apply();
		}

		// DEBUG
		//if( hasTraComments() ) getTraComments()->setFormat( "CLUSTER-REP: Cluster-ene is: %8.3lf. CRMS: %8.3lf" )
		//	(loopSet.posCachePointer[i].first / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na))
		//	((*(loopSet.posCachePointer[i].second)).calcCRMS(PickHeavyAtoms()));				
		//getWSpace().outtra.append();
		// END DEBUG

		m_Restraint->activate();
		m_Restraint->forceSetup();
		getWSpace().Step = 0;
		m_Min->run();
		m_Restraint->deactivate();
		getWSpace().Step = 0;
		m_Min->run();

		double ene = getWSpace().ene.epot / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na);

		// DEBUG
		//if( hasTraComments() ) getTraComments()->setFormat( "GBSA-Minimised CLUSTER-REP: Minimised-ene is: %8.3lf. CRMS: %8.3lf" )
		//	(ene)
		//	((*(loopSet.posCachePointer[i].second)).calcCRMS(PickHeavyAtoms()));				
		//getWSpace().outtra.append();
		// END DEBUG

		getWSpace().Step = outerStep + 1;  // important :-D	

		return ene;
	}

	void ArcusRefine_ToTra::setup( LoopSet& loopSet )
	{
		StringBuilder filename( m_FileStem );
		filename.appendFormat( "_%d_%d.tra" )(loopSet.range.getStartResIndex())(loopSet.range.getNRes());
		m_Tra = new IO::OutTra_BTF( filename.toString(), getWSpace() );
		m_Tra->create();
	}

	double ArcusRefine_ToTra::refine( LoopSet& loopSet, size_t _ClusterRepIndex )
	{
		if( OutputLevel >= Verbosity::Normal )
		{
			Printf("Invoking ArcusRefine_ToTra:\n");
		}

		m_WSpace->Step = _ClusterRepIndex;

		size_t i = loopSet.clusterReps[_ClusterRepIndex];
		loopSet.posCachePointer[i].second->revert();

		m_FF->calcForces(); // Ensure the tra contains the correct energy information

		double cRMS = getWSpace().calcCRMS_HeavyAtom(); // and the correct cRMS value
		double epot = getWSpace().ene.epot;

		if( hasTraComments() ) 
		{
			getTraComments()->setFormat(
				"Whole protein heavy-atom cRMS: %6.3lf, Loop mainchain cRMS: %6.3lf, Loop heavy-atom cRMS: %6.3lf, Epot: %8.3lf")
				(cRMS)
				(loopSet.posCachePointer[i].second->calcCRMS( PickCoreBackbone() ) )
				(loopSet.posCachePointer[i].second->calcCRMS( PickHeavyAtoms() ) )
				( epot / Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na );
		}
		m_Tra->append();	

		return epot;
	}
}

