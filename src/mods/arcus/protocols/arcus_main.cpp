#include "global.h"
#include "tools/vector.h"
#include "tools/statclock.h"
#include "library/angleset.h"
#include "library/rotamerlib.h"
#include "workspace/workspace.h" 
#include "workspace/neighbourlist.h"
#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "forcefields/segrejoin.h"
#include "protocols/torsionalminimisation.h"
#include "manipulators/rotamer_scwrl.h"
#include "filters/loopcadistfilter.h"
#include "tools/vector.h"
#include "arcus.h"

#define DIST_FILT // enable distance filtering code
#define SEG_FILT // enable segment-join filtering code

namespace Protocol
{
	Arcus::Arcus( Physics::Forcefield& _ffs, Physics::Forcefield& _ff, const Library::AngleSet& _as, const Library::RotamerLibrary& _rotLib )
		: ArcusBase( _ffs, _ff, _as ), 
		m_RotLib(&_rotLib), 
		m_StaticGrid4A( _ffs.getWSpace(), PickNothing(), 0.0 ), // NULL, shite state for init
		m_StaticGrid6_5A( _ffs.getWSpace(), PickNothing(), 0.0 ) // NULL, shite state for init
	{
		setDefaultRefiner();
		ProduceNJoinedBranches = SIZE_T_FAIL; // max value, equivelent to disabled.
		m_SplitBuilder = true; // pretty much fundamental to this ArcusBase derived implementation
	}

	void Arcus::calibrateBuilders()
	{
		for( size_t i = 0; i < m_ConfBuilder.size(); i++ )
		{
			m_ConfBuilder[i].PropensityWeighting = true;
			m_ConfBuilder[i].NoiseSigma = Maths::DegToRad(15.0);
			m_ConfBuilder[i].setRandCount(SIZE_T_FAIL, false); // this sets infinity, we usually dont need this many ;-)
			m_ConfBuilder[i].setBuildMode(Manipulator::ConfBuilderBase::Backbone);
		}
	}

	void Arcus::configLevel1Filters( size_t i )
	{
		// Make the relevent filters for each segment
		m_Filters.push_back( FilterContainer() );

		// 1) Self clash filter
		ClashFilter* filterSelfClash = new ClashFilter(); 
		filterSelfClash->setOLapFac(0.65); // Mild clashes are allowed
		PickCoreBackbone bbOnly;
		Pick_AND pickSelf(PickResidueRange(m_ConfBuilderSeg[i]),bbOnly);
		filterSelfClash->setForPicker(pickSelf);
		filterSelfClash->setAgainstPicker(pickSelf); // Against all non-moving atoms
		filterSelfClash->setMolecule( getWSpace() );
		m_Filters[i].addWithOwnership( filterSelfClash );

		// 2) SegmentDistanceFilter: Prevent over-extension
#ifdef SEG_FILT
		if( SegFanFilename.size() > 0 )
		{
			ASSERT( IO::fileExists( SegFanFilename ), IOException, "Cannot find seg-dat file");
			SegmentDistanceFilter* distFilter = new SegmentDistanceFilter(SegFanFilename);
			distFilter->initialise( m_Regions[i/2], i % 2 != 0 ); // i % 2 != 0 -- detects N vs C terminal builder
			m_Filters[i].addWithOwnership( distFilter );
		}
#endif

		// 3) Protein body distance filter
#ifdef DIST_FILT
		size_t NOffset = 0; 
		size_t COffset = 0;
		if( i % 2 == 0 )
		{
			NOffset++; // Efficiency - don't bother with the residue next to the anchor...
		}
		else
		{
			COffset++; // Efficiency - don't bother with the residue next to the anchor...
		}
		LoopCADistFilter* distFromProtinBody = new LoopCADistFilter();
		distFromProtinBody->setTo( m_StaticGrid6_5A, getWSpace(), PickResidueRange(m_ConfBuilderSeg[i]), NOffset, COffset );
		m_Filters[i].addWithOwnership( distFromProtinBody );
#endif

		// 4) Surface clash filter
		ClashGridFilter* filter = new ClashGridFilter(); 
		filter->setOLapFac(0.65); // Mild clashes are allowed
		filter->setForPicker( pickSelf );
		filter->setAgainstPicker( Pick_AND(Pick_NOT(m_DynamicAtomsPicker),PickHeavyAtoms()) ); // Against all non-moving atoms
		filter->setMolecule( getWSpace() );
		m_Filters[i].addWithOwnership( filter );
	}

	void Arcus::generateValidConformers()
	{
		ASSERT( m_RotLib != NULL, NullInternalException, "Internal NULL for RotLib in Arcus");
		ASSERT( m_ConfBuilderSeg.size() == m_Regions.size() * 2, CodeException, "Arcus needs this to be true");

		// We are using CoreBackbone perturbation below, meaning that the H and HA atoms and the sidechain dont move
		// Once that has passes, the RotamerApplication will apply the sidechain, but the backbone Hydrogens still
		// require placement. This is performed by this class, which forces rebuild for all atoms that match the pattern
		// isBackbone() && isHydrogen(), **regardless** of isRebuildRequired();
		MainchainHydrogenRebuilder bbHydrogenBuilder;

		// Set the conformer builders to do what we want for this derived class
		calibrateBuilders();

		// ---------------------
		//  Stage 1 - Main Loop
		// ---------------------

		const size_t branchRepeatAlloc = 50;
	
		m_BuiltSections.resize(m_Regions.size());
		
		for( size_t q = 0; q < m_Regions.size(); q++ )
		{
			const size_t IndexerN = q*2;
			const size_t IndexerC = q*2+1;

			// -----------------------
			//  Stage 1A - EQUIPEMENT
			// -----------------------

			size_t currentBranchAlloc = 0;

			std::vector<PosStore> NBranches;
			std::vector<PosStore> CBranches;
			// Reserve enough space to avoid constant reallocation.
			size_t branchReserve = ProduceNJoinedBranches * 100; // 100 -> 1% of branch-pairs can be joined - ish...
			NBranches.reserve( branchReserve );
			CBranches.reserve( branchReserve );

			// 4) N-terminal Surface clash filter, now WITH sidechains
			ClashGridFilter filterN; 
			filterN.GridCutoff = 3.5; // Atoms within 3.5 of a grid are included in the AGAINST list
			filterN.setOLapFac(0.65); // Modrate clashes are allowed
			filterN.setForPicker( Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerN]),PickHeavyAtoms()) );
			filterN.setAgainstPicker( Pick_AND(Pick_NOT(m_DynamicAtomsPicker),PickHeavyAtoms()) ); // Against all non-moving atoms
			filterN.setMolecule( getWSpace() );

			// 4) C-terminal Surface clash filter, now WITH sidechains
			ClashGridFilter filterC; 
			filterC.GridCutoff = 3.5; // Atoms within 3.5 of a grid are included in the AGAINST list
			filterC.setOLapFac(0.65); // Modrate clashes are allowed
			filterC.setForPicker( Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerC]),PickHeavyAtoms()) );
			filterC.setAgainstPicker( Pick_AND(Pick_NOT(m_DynamicAtomsPicker),PickHeavyAtoms()) ); // Against all non-moving atoms
			filterC.setMolecule( getWSpace() );

			Manipulator::RotamerApplicator_SCWRL rotN( getWSpace(), *m_RotLib, 
				PickResidueList( m_ConfBuilderSeg[IndexerN] )
				);
			rotN.SupressStaticGridRefreshOnApply = true;
			rotN.OutputLevel = Verbosity::Silent;
			rotN.overrideStaticGrid(m_StaticGrid4A); // rotamers can only now clash with the rigid-body

			Manipulator::RotamerApplicator_SCWRL rotC( getWSpace(), *m_RotLib, 
				PickResidueList( m_ConfBuilderSeg[IndexerC] )
				);
			rotC.SupressStaticGridRefreshOnApply = true;
			rotC.OutputLevel = Verbosity::Silent;
			rotC.overrideStaticGrid(m_StaticGrid4A); // rotamers can only now clash with the rigid-body

			double bestcrmsN = DBL_MAX;
			double bestcrmsC = DBL_MAX;

			// Make a currentPos so that internal setup happends only once
			// We will add copies of this to the conatiner. Copies are cheap..
			PosStore currentBranchN( getWSpace(), m_ConfBuilderSeg[IndexerN] );
			PosStore currentBranchC( getWSpace(), m_ConfBuilderSeg[IndexerC] );

			int cyclesN = 0;
			int cyclesPassedN = 0;
			int cyclesC = 0;
			int cyclesPassedC = 0;

			// -----------------------
			//  Stage 1B - EQUIPEMENT
			// -----------------------

			long cycles = 0;
			long cyclesPassed = 0;

			std::vector<double> bbcrms;
			std::vector<double> aacrms;

			// Make a currentPos so that internal setup happens only once
			// We will add copies of this to the conatiner. Copies are cheap..
			PosStore curentPos( getWSpace(), m_Regions[q] );
			// Reserve enough space to avoid constant reallocation. 1 in 1000 is optomistic...

			LoopSet& loopSet = m_BuiltSections[q];
			loopSet.range = PickResidueRange( m_Regions[q] );
			loopSet.posCache.reserve( ProduceNJoinedBranches );

			// 1) Join Filter
			OmegaGroupFilter joinFilter;
			joinFilter.AssessSecondaryDistances = false; // At the moment we are only interested in the ball-park
			joinFilter.AssessTorsion = false; // At the moment we are only interested in the ball-park
			joinFilter.setTo( getWSpace(), m_Regions[q].getBreakResIndex(), 1.5 ); // 1.5Angstrom join diustance for all atom-pairs

			// 2) Clash-between-NC-pair filter
			ClashFilter filterPairwiseClash; 
			filterPairwiseClash.setOLapFac(0.65); // Mild clashes are allowed
			filterPairwiseClash.setForPicker(Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerN]),PickCoreBackbone()));
			filterPairwiseClash.setAgainstPicker(Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerC]),PickCoreBackbone()));
			filterPairwiseClash.setMolecule( getWSpace() );

			Manipulator::RotamerApplicator_SCWRL rotBoth( getWSpace(), *m_RotLib, 
				PickResidueList( m_Regions[q] )
				);
			rotBoth.SupressStaticGridRefreshOnApply = true; // we know the core isn't moving. Disable change-detection.
			rotBoth.OutputLevel = Verbosity::Silent;
			rotBoth.overrideStaticGrid(m_StaticGrid4A); // rotamers can only now clash with the rigid-body

			// 3) Self clash filter, HEAVY atoms now sidechains are packed
			ClashFilter filterPairwiseClash2; 
			filterPairwiseClash2.setOLapFac(0.65); // Mild clashes are allowed
			filterPairwiseClash2.setForPicker(Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerN]),PickHeavyAtoms()));
			filterPairwiseClash2.setAgainstPicker(Pick_AND(PickResidueRange(m_ConfBuilderSeg[IndexerC]),PickHeavyAtoms()));
			filterPairwiseClash2.setMolecule( getWSpace() );

			// 4) Surface clash filter, now WITH sidechains
			ClashGridFilter finalSurfaceFilter; 
			finalSurfaceFilter.GridCutoff = 3.5; // Atoms within 3.5 of a grid are included in the AGAINST list
			finalSurfaceFilter.setOLapFac(0.65); // Modrate clashes are allowed
			finalSurfaceFilter.setForPicker( Pick_AND(PickResidueRange(m_Regions[q]),PickHeavyAtoms()) );
			finalSurfaceFilter.setAgainstPicker( Pick_AND(Pick_NOT(m_DynamicAtomsPicker),PickHeavyAtoms()) ); // Against all non-moving atoms
			finalSurfaceFilter.setMolecule( getWSpace() );

			// -----------------------
			//  Stage 1 - Uber-LOOP
			// -----------------------
			do
			{
				// ---------------------
				//  Stage 1a - half-gen
				//  N-Terminus
				// ---------------------

				if( OutputLevel >= Verbosity::Normal )
					Printf("N-terminal enumeration for N-conformers %d to %d\n")
					(currentBranchAlloc)(currentBranchAlloc+branchRepeatAlloc);

				size_t innerPassed = 0;
				while( innerPassed < branchRepeatAlloc )
				{
					m_ConfBuilder[IndexerN].next();
					cyclesN++;
					if( m_Filters[IndexerN].passes() )
					{		
						// This call actually turns the BB-only model into an all-atom model - neat hu, :-D
						rotN.apply();
						//if( hasTraComments() ) getTraComments()->setTo( "Packed rot..." );
						//getWSpace().outtra.append();

						// Ensure that the H and HA atoms are in sensible places too!!
						ASSERT( bbHydrogenBuilder.invokeBuild( getWSpace(), Verbosity::Silent, PickResidueList( m_ConfBuilderSeg[IndexerN] ), m_ConfBuilderSeg[IndexerN] ),
							ProcedureException, "MainchainHydrogenRebuilder failed on internal call" );
						//if( hasTraComments() ) getTraComments()->setTo( "Built BB hydrogens..." );
						//getWSpace().outtra.append();

						if( filterN.passes() ) // test the sidechains are packed well enough
						{		
							innerPassed++;
							cyclesPassedN++;
							currentBranchN.store();
							NBranches.push_back(currentBranchN);
							double crmsN = m_ConfBuilderSeg[IndexerN].cCRMS_BB();
							if( bestcrmsN > crmsN ) bestcrmsN = crmsN;
							//if( hasTraComments() )getTraComments()->setFormat( "Passed :-D: backbone crms:%8.3lf, best_crms:%8.3lf\n" )(crms)(bestcrms);
							//getWSpace().outtra.append();
						}
						else
						{
							//std::cout << filter.reason();
							//if( hasTraComments() ) getTraComments()->setTo( filter.reason() );
							//getWSpace().outtra.append();
						}
					}
					else
					{
						//std::cout << m_Filters[j].reason();
						//if( hasTraComments() ) getTraComments()->setTo( m_Filters[j].reason() );
						//getWSpace().outtra.append();
					}
				}

				if( OutputLevel >= Verbosity::Normal )
					Printf("So far --> Best bb-crms (N-builder %d) was %5.3lf, pass %6.3lf%%, (%d / %d attempted)\n")
						(IndexerN)(bestcrmsN)( 100.0 * ((double)cyclesPassedN/(double)cyclesN) )(cyclesPassedN)(cyclesN);

				// ---------------------
				//  Stage 1a - half gen
				//  C-Terminus
				// ---------------------

				if( OutputLevel >= Verbosity::Normal )
					Printf("C-terminal enumeration for C-conformers %d to %d\n")
					(currentBranchAlloc)(currentBranchAlloc+branchRepeatAlloc);

				innerPassed = 0;
				while( innerPassed < branchRepeatAlloc )
				{
					m_ConfBuilder[IndexerC].next();
					cyclesC++;
					if( m_Filters[IndexerC].passes() )
					{		
						// This call actually turns the BB-only model into an all-atom model - neat hu, :-D
						rotC.apply();
						//if( hasTraComments() ) getTraComments()->setTo( "Packed rot..." );
						//getWSpace().outtra.append();

						// Ensure that the H and HA atoms are in sensible places too!!
						ASSERT( bbHydrogenBuilder.invokeBuild( getWSpace(), Verbosity::Silent, PickResidueList( m_ConfBuilderSeg[IndexerC] ), m_ConfBuilderSeg[IndexerC] ),
							ProcedureException, "MainchainHydrogenRebuilder failed on internal call" );
						//if( hasTraComments() ) getTraComments()->setTo( "Built BB hydrogens..." );
						//getWSpace().outtra.append();

						if( filterC.passes() ) // test the sidechains are packed well enough
						{		
							innerPassed++;
							cyclesPassedC++;
							currentBranchC.store();
							CBranches.push_back(currentBranchC);
							double crmsC = m_ConfBuilderSeg[IndexerC].cCRMS_BB();
							if( bestcrmsC > crmsC ) bestcrmsC = crmsC;
							//if( hasTraComments() )getTraComments()->setFormat( "Passed :-D: backbone crms:%8.3lf, best_crms:%8.3lf\n" )(crms)(bestcrms);
							//getWSpace().outtra.append();
						}
						else
						{
							//std::cout << filter.reason();
							//if( hasTraComments() ) getTraComments()->setTo( filter.reason() );
							//getWSpace().outtra.append();
						}
					}
					else
					{
						//std::cout << m_Filters[j].reason();
						//if( hasTraComments() ) getTraComments()->setTo( m_Filters[j].reason() );
						//getWSpace().outtra.append();
					}
				}

				if( OutputLevel >= Verbosity::Normal )
					Printf("So far --> Best bb-crms (C-builder %d) was %5.3lf, pass %6.3lf%%, (%d / %d attempted)\n")
						(IndexerC)(bestcrmsC)( 100.0 * ((double)cyclesPassedC/(double)cyclesC) )(cyclesPassedC)(cyclesC);

				// -----------------------
				// Stage 1b - half pairing
				// -----------------------

				if( OutputLevel >= Verbosity::Normal )
					Printf("Pairing N and C branches...\n");

				ASSERT( NBranches.size() - currentBranchAlloc  ==  branchRepeatAlloc, CodeException, "Assumption Failure" );
				ASSERT( NBranches.size() == CBranches.size(), CodeException, "Assumption Failure" );

				for( size_t k = 0; k < NBranches.size(); k++ )
				{
					NBranches[k].revert();

#ifdef DIST_FILT
					// REQUIRED FOR **DEBUG** only
					//LoopCADistFilter distFromProtinBody;
					//distFromProtinBody.setTo( m_StaticGrid6_5A, getWSpace(), PickResidueRange(m_ConfBuilderSeg[IndexerN]), 1, 0 );
					//ASSERT( distFromProtinBody.passes(), CodeException, "Internal validation failure" );
					// END DEBUG only
#endif

					// Decide which C-branches to use - we dont want ones compared in previous sweeps!
					// They have already been assessed as pairs!!
					size_t CFrom = ( k >= currentBranchAlloc ) ? 0 : currentBranchAlloc;				
					size_t CTo = CBranches.size();
					for( size_t m = CFrom; m < CTo; m++ )
					{
						cycles++;
						CBranches[m].revert();

#ifdef DIST_FILT
						// REQUIRED FOR **DEBUG** only
						//LoopCADistFilter distFromProtinBody;
						//distFromProtinBody.setTo( m_StaticGrid6_5A, getWSpace(), PickResidueRange(m_ConfBuilderSeg[IndexerC]), 0, 1 );
						//ASSERT( distFromProtinBody.passes(), CodeException, "Internal validation failure" );
						// END DEBUG only
#endif

						if( joinFilter.passes() )
						{
							if( filterPairwiseClash.passes() )
							{ 
								// Apply rotamers, this time to the whole conformer
								rotBoth.apply();
								if( filterPairwiseClash2.passes() )
								{
									if( finalSurfaceFilter.passes() )
									{
										cyclesPassed++;

										// The structure has passed all initial tests.
										// Store it, and add to the current posCache container
										curentPos.store();
										loopSet.posCache.push_back(curentPos);
 										bbcrms.push_back(curentPos.calcCRMS( PickCoreBackbone() ) );
										aacrms.push_back( curentPos.calcCRMS( PickHeavyAtoms() ) );

										if( m_BuiltSections[q].posCache.size() == ProduceNJoinedBranches )
											goto COOKED;
										 
										//if( hasTraComments() ) getTraComments()->setFormat( "Niice: %8.3lf, best was %8.3lf\n" )(crms)(bestcrms);
										//getWSpace().outtra.append();

#ifdef DIST_FILT
										// REQUIRED FOR **DEBUG** only
										//LoopCADistFilter distFromProtinBody;
										//distFromProtinBody.setTo( m_StaticGrid6_5A, getWSpace(), m_Regions[q], 1, 1 );
										//if( !distFromProtinBody.passes() )
										//{
										//	std::string reason = distFromProtinBody.reason();
										//	std::cout << reason;
										//	if( hasTraComments() ) getTraComments()->setTo( reason );
										//	getWSpace().outtra.append();
										//	THROW( CodeException, "Internal validation failure!");
										//}
										// END DEBUG
#endif
									}
									else
									{
										//if( hasTraComments() ) getTraComments()->setTo( finalSurfaceFilter.reason() );
										//getWSpace().outtra.append();
									}
								}
								//else
								//{
								//	if( hasTraComments() ) getTraComments()->setTo( filterPairwiseClash2.reason() );
								//	getWSpace().outtra.append();
								//}
							}
							//else
							//{
							//	if( hasTraComments() ) getTraComments()->setTo( filterPairwiseClash.reason() );
							//	getWSpace().outtra.append();
							//}
						}
						//else
						//{
						//	if( hasTraComments() ) getTraComments()->setTo( joinFilter.reason() );
						//	getWSpace().outtra.append();
						//}
					}
				}

				if( OutputLevel >= Verbosity::Normal )
					Printf("%d pairs so far!\n")(m_BuiltSections[q].posCache.size() );

				// -----------------------
				// Stage-1 loop end
				// Rather important!!
				// We dont want to compare the same branch-pair twice!!!
				// -----------------------
				currentBranchAlloc += branchRepeatAlloc;				
			}
			while( m_BuiltSections[q].posCache.size() < ProduceNJoinedBranches );

COOKED:
			// --------------------------------------
			//  Final Stage-1 Reporting, per-section
			// --------------------------------------

			Printf("Region %d report:\n")(q);
			Printf("Best bb-crms (N-builder %d) was %5.3lf, pass %6.3lf%%, (%d / %d attempted)\n")
				(IndexerN)(bestcrmsN)( 100.0 * ((double)cyclesPassedN/(double)cyclesN) )(cyclesPassedN)(cyclesN);
			Printf("Best bb-crms (C-builder %d) was %5.3lf, pass %6.3lf%%, (%d / %d attempted)\n")
				(IndexerC)(bestcrmsC)( 100.0 * ((double)cyclesPassedC/(double)cyclesC) )(cyclesPassedC)(cyclesC);

			double bestBB = Min(bbcrms);
			double bestAA = Min(aacrms);

			Printf("Best rejoin-pair (%d) backbone crms:%5.3lf heavy-atom crms:%5.3lf, pass %6.3lf%% (%ld / %ld attempted)\n")
				(q)
				(bestBB)
				(bestAA)
				( 100.0 * ((double)cyclesPassed/(double)cycles) )
				(cyclesPassed)
				(cycles);

			// *************************************
			// UBER DEBUG !!!!
			//StringBuilder sb(nameMe);
			//sb.appendFormat(",%d")(ProduceNJoinedBranches);
			//if( bbcrms.size() != 0 )
			//{
			//	sb.appendFormat(",%8.3lf,%8.3lf,%8.3lf,%8.3lf,%8.3lf,%8.3lf,%6.3lf,%ld,%ld\n")
			//		( Average(bbcrms) )( Average(aacrms) )
			//		( bestBB )( bestAA )
			//		( Max(bbcrms) )( Max(aacrms) )
			//		( 100.0 * ((double)cyclesPassed/(double)cycles) )
			//		(cyclesPassed)
			//		(cycles);
			//}
			//else
			//{
			//	sb.append(",Nonefound\n");
			//	Printf(sb.toString());
			//}
			//std::string writeTo = nameMe + ".csv";
			//std::ofstream write( writeTo.c_str(), std::ifstream::out );
			//ASSERT( write.is_open(), CodeException, "Anusol");
			//write << sb;
			//write.close();
			// UBER DEBUG !!!!
			// *************************************
		}

		// DEBUG-only reporter for-loop ....
		//for( size_t j = 0; j < m_Regions.size(); j++ )
		//{
		//	LoopSet& loopSet = m_BuiltSections[j];
		//	std::vector< PosStore >& joinCapablePairs = loopSet.posCache;
		//	for( size_t i = 0; i < joinCapablePairs.size(); i++ )
		//	{
		//		joinCapablePairs[i].revert();
		//		if( hasTraComments() ) getTraComments()->setFormat( "Join-Potential structure heavy-atom CRMS: %8.3lf" )
		//			(joinCapablePairs[i].calcCRMS(PickHeavyAtoms()));
		//		getWSpace().outtra.append();
		//	}
		//}

		return;
	}

	void Arcus::rejoinConformers()
	{
		// ----------------------------
		//  Stage 2a - half re-joining
		// ----------------------------

		OmegaGroupFilter joinFilter;
		joinFilter.AssessTorsion = false;

		AngleSetFilter anglesetFilter;	

		LoopCADistFilter distFromProtinBody;

		for( size_t j = 0; j < m_Regions.size(); j++ )
		{
			joinFilter.setTo( getWSpace(), m_Regions[j].getBreakResIndex() );
			anglesetFilter.setTo( getWSpace(), getAngleSet(), PickResidueList(m_Regions[j]) );
			distFromProtinBody.setTo( m_StaticGrid6_5A, getWSpace(), m_Regions[j], 1, 1 );

			LoopSet& loopSet = m_BuiltSections[j];
			std::vector< PosStore >& acceptedPairsSet = loopSet.posCache;
			size_t nStage1Find = acceptedPairsSet.size();
			if( nStage1Find == 0 ) 
				throw ProcedureException("Stage-2 has nothing to work with mein-herring!!");
			std::vector< std::pair<double,PosStore*> >& acceptedPairsEne = loopSet.posCachePointer;
			acceptedPairsEne.reserve( acceptedPairsSet.size() );

			Physics::PeptideGroupRejoinForce rejoinForce( m_Regions[j] );
			rejoinForce.Initialise( 50.0 ); // 50 does the trick...
			ffs->add( rejoinForce );

			PickAtomRanges atRange;
			setPickedRegions_Core( getWSpace(), atRange, m_Regions[j] );
			Protocol::TorsionalMinimisation torMin( *ffs, atRange );
			torMin.Steps = 1000;
			torMin.UpdateNList = UpdateNList; // propogate the parameter to the child-protocol
			torMin.SlopeCutoff = 0.1 * Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na;
			torMin.InitialCapFactor = 0.05;
			torMin.OutputLevel = Verbosity::Silent;
			torMin.UpdateMon = 0;
			//torMin.UpdateScr = 50;
			//torMin.UpdateTra = 10;

			for( size_t k = 0; k < acceptedPairsSet.size(); k++ )
			{
				acceptedPairsSet[k].revert();

				// Stop the backbone buggering off...
				rejoinForce.enableHarmonicCARestraints(true);

#ifdef DIST_FILT
				// This should have been checked prior to minimisation ...
				if( !distFromProtinBody.passes() )
				{
					std::string reason = distFromProtinBody.reason();
					std::cout << reason;
					if( hasTraComments() ) getTraComments()->setTo( reason );
					getWSpace().outtra.append();
					THROW( CodeException, "Internal validation failure!");
				}
#endif

				torMin.run();

				if( !joinFilter.passes() )
				{
					if( OutputLevel >= Verbosity::Loud && hasTraComments() ) 
					{
						getTraComments()->setTo( "Steric tor-min failed to re-join\n" );
						getTraComments()->toScreen();
						getWSpace().outtra.append();
					}
					continue;
				}

				if( !distFromProtinBody.passes() )
				{
					if( OutputLevel >= Verbosity::Loud && hasTraComments() ) 
					{
						getTraComments()->setTo( "Steric tor-min failed ensure protein-surface-contact\n" );
						getTraComments()->toScreen();
						getWSpace().outtra.append();
					}
					continue;
				}

				if( !anglesetFilter.passes() )
				{
					if( OutputLevel >= Verbosity::Loud && hasTraComments() ) 
					{
						getTraComments()->setTo( "Steric tor-min failed to retain valid torsions\n" );
						getTraComments()->toScreen();
						getWSpace().outtra.append();
					}
					continue;
				}

				acceptedPairsSet[k].store();
				acceptedPairsEne.push_back( std::pair<double,PosStore*>( getWSpace().ene.epot, &acceptedPairsSet[k] ) );
				if( hasTraComments() ) 
				{
					getTraComments()->setFormat( "Cool structure :-D Ene: %8.3lf. Heavy-Atom CRMS: %8.3lf\n" )
						(getWSpace().ene.epot / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na) )
						(acceptedPairsSet[k].calcCRMS(PickHeavyAtoms()));
					getTraComments()->toScreen();
				}
				getWSpace().outtra.append();

				continue;
			}

			ffs->pop_back();
		}

		// DEBUG-only reporter for-loop ....
		//for( size_t j = 0; j < m_Regions.size(); j++ )
		//{
		//	LoopSet& loopSet = m_BuiltSections[j];
		//	std::vector< std::pair<double,PosStore*> >& acceptedPairsEne = loopSet.posCachePointer;
		//	for( size_t i = 0; i < acceptedPairsEne.size(); i++ )
		//	{
		//		if( hasTraComments() ) getTraComments()->setFormat( "Cluster-ene is: %8.3lf. Heavy-Atom CRMS: %8.3lf" )
		//			(acceptedPairsEne[i].first / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na))
		//			((*(acceptedPairsEne[i].second)).calcCRMS(PickHeavyAtoms()));
		//		acceptedPairsEne[i].second->revert();
		//		getWSpace().outtra.append();
		//	}
		//}

		// ----------------------------------------------
		//  Stage 2b - clustering
		//  or more specifically similarity filtering...
		// ----------------------------------------------

		const double RMS_Cutoff = 0.25; // may be sensible ?

		for( size_t j = 0; j < m_Regions.size(); j++ )
		{
			// Proxies
			LoopSet& loopSet = m_BuiltSections[j];
			std::vector< PosStore >& acceptedPairsSet = loopSet.posCache;
			std::vector< std::pair<double,PosStore*> >& acceptedPairsEne = loopSet.posCachePointer;
			std::vector<size_t>& clusterReps = loopSet.clusterReps;

			// Critical for the algorithm below to work!
			std::sort( acceptedPairsEne.begin(), acceptedPairsEne.end() );

			// We want this info
			std::vector<bool> used( acceptedPairsEne.size(), false );
			size_t sweepStart = 0;
			size_t currentRep;
			while(true)
			{
				currentRep = SIZE_T_FAIL;
				for( size_t i = sweepStart; i < used.size(); i++ )
				{
					if( used[i] ) 
						continue;
					used[i] = true;
					currentRep = i;
					clusterReps.push_back(i);
					break;
				}
				if( currentRep == SIZE_T_FAIL )
					break; // we have traversed the list. All allocations have been made
				sweepStart = currentRep + 1;
				PosStore& rep = *acceptedPairsEne[currentRep].second;

				for( size_t k = sweepStart; k < used.size(); k++ )
				{
					if( used[k] ) 
						continue;

					PosStore& posK = *acceptedPairsEne[k].second;

					double crms = rep.calcCRMSOfStoreTo( posK, PickHeavyAtoms() );
					if( crms < RMS_Cutoff )
					{
						used[k] = true;
						// If we are making clusters, add it here...
					}
				}
			}
		}

		// DEBUG info...
		//for( size_t j = 0; j < m_Regions.size(); j++ )
		//{
		//	// Proxies
		//	LoopSet& loopSet = m_BuiltSections[j];
		//	std::vector<size_t>& clusterReps = loopSet.clusterReps;
		//	std::vector< std::pair<double,PosStore*> >& acceptedPairsEne = loopSet.posCachePointer;
		//	for( size_t k = 0; k < clusterReps.size(); k++ )
		//	{
		//		size_t i = clusterReps[k];
		//		if( hasTraComments() ) getTraComments()->setFormat( "CLUSTER-REP: Cluster-ene is: %8.3lf. Heavy-Atom CRMS: %8.3lf" )
		//			(acceptedPairsEne[i].first / (Physics::PhysicsConst::kcal2J / Physics::PhysicsConst::Na))
		//			((*(acceptedPairsEne[i].second)).calcCRMS(PickHeavyAtoms()));
		//		acceptedPairsEne[i].second->revert();
		//		getWSpace().outtra.append();
		//	}
		//}

		return;
	}

	int Arcus::initialise()
	{
		// What are we actually moving?
		setPickedRegions( getWSpace(), m_DynamicAtomsPicker, m_Regions );
		m_StaticGrid4A = ProximityGrid(getWSpace(),Pick_NOT(m_DynamicAtomsPicker),4.0);
		m_StaticGrid6_5A = ProximityGrid(getWSpace(),Pick_NOT(m_DynamicAtomsPicker),6.50);

		return ArcusBase::initialise();
	}

	void Arcus::AssertCA6_32Filter()
	{
		m_StaticGrid6_5A = ProximityGrid(getWSpace(),Pick_NOT(m_DynamicAtomsPicker),6.50);
		m_StaticGrid6_5A.refresh();
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			LoopCADistFilter checkThis;
			checkThis.setTo( m_StaticGrid6_5A, getWSpace(), m_Regions[i], 1, 1 );
			getWSpace().outtra.append();
			bool check = checkThis.passes();
			if( !check ) std::cout << checkThis.reason();
			ASSERT( check, CodeException, "AssertCA6_32Filter fail!");
		}
	}

	void Arcus::AssertSegDistFilter()
	{
		ASSERT( SegFanFilename.size() > 0, NullInternalException, "You can't test the SegDist filter using 'AssertSegDistFilter()' if 'SegFanFilename' is undefined...");
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			SegmentDistanceFilter filter( SegFanFilename );
			filter.initialise( m_Regions[i], false ); // N-terminal, forware-orientation.
			ASSERT( filter.passes(), CodeException, filter.reason() );			
			filter.initialise( m_Regions[i], true ); // C-terminal, reversed-orientation.
			ASSERT( filter.passes(), CodeException, filter.reason() );	
		}
	}

	void Arcus::setDefaultRefiner()
	{
		m_Stage3Refiner = &m_DefaultStage3;
	}

	void Arcus::setRefiner( ArcusRefineBase& _refiner )
	{
		m_Stage3Refiner = &_refiner;
	}

	int Arcus::runcore()
	{
		// Re-Initialise
		m_BuiltSections.clear();
		m_StaticGrid4A.refresh();
		m_StaticGrid6_5A.refresh();

		StatClock clockMe;
		clockMe.setName("Stage-1");
		clockMe.Begin();

		// ***********************************
		// Fill the m_BuiltSections array with stuff. 
		// You never know, it might even be decent!
		if( OutputLevel >= Verbosity::Normal ) 
			Printf("Beginning conformer selection\n");		
		generateValidConformers(); // Stage-1
		// ***********************************

		clockMe.End();

		if( OutputLevel >= Verbosity::Normal ) 
		{
			Printf("Done!\n");
			clockMe.ReportSeconds();
			Printf("Valid conformer pairs found:\n");
			for( size_t i = 0; i < m_BuiltSections.size(); i++ )
			{
				Printf("  Section %d (from %d to %d) has %d viable conformations\n")(i)
					(m_Regions[i].getStartResIndex())
					(m_Regions[i].getEndResIndex())
					(m_BuiltSections[i].posCache.size());
			}
			Printf("Beginning conformer rejoin\n");
		}

		clockMe.Reset();
		clockMe.setName("Stage-2");
		clockMe.Begin();

		// ****************************
		rejoinConformers(); // Stage-2
		// ****************************

		clockMe.End();

		if( OutputLevel >= Verbosity::Normal ) 
		{
			Printf("Done!\n");
			clockMe.ReportSeconds();
			Printf("Valid conformer pairs found:\n");
			for( size_t i = 0; i < m_BuiltSections.size(); i++ )
			{
				Printf("  Section %d (from %d to %d) has %d clustered conformations\n")(i)
					(m_Regions[i].getStartResIndex())
					(m_Regions[i].getEndResIndex())
					(m_BuiltSections[i].clusterReps.size());
			}
			Printf("Beginning conformer refinement\n");
		}

		clockMe.Reset();
		clockMe.setName("Stage-3");
		clockMe.Begin();

		// ***********************************
		m_Stage3Refiner->enableComments( getTraComments() );
		m_Stage3Refiner->initialise( getWSpace(), *ff );
		m_Stage3Refiner->refine( m_BuiltSections ); // Stage-3 (bridge-pattern)
		// ***********************************

		clockMe.End();

		if( OutputLevel >= Verbosity::Normal ) 
		{
			Printf("Done!\n");
			clockMe.ReportSeconds();
		}

		return 1;
	}
}

