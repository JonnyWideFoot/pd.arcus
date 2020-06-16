#include "global.h"

#include <sstream> // String stream
#include "tools/vector.h"
#include "library/backbonetorsions.h"
#include "library/angleset.h"
#include "segtorsions.h"
#include "system/segbuilder.h"

#include "forcefields/breakablebonded.h"
#include "protocols/minimise.h"

// Self Header Include
#include "conformerbuilderbase.h"

using namespace std;
using namespace Maths;
using namespace Library; ///< so that we can use the 'conformer_type' definition

namespace Manipulator
{

	// -----------------
	//  ConfBuilderBase
	// -----------------

	ConfBuilderBase::ConfBuilderBase( const Library::AngleSet& _angleSet, SegmentDef& _segdef ):
		SegmentDefUser(_segdef),
		SegCoreBBOnlyTC(_segdef),
		SegWholeBackboneTC(_segdef),
		m_AngleSet( &_angleSet ),
		OutputLevel( Verbosity::Normal ),
		m_IdealisedMinimised( false ),
		m_BuildMode(AllAtom), // The safest should be the default
		NoiseSigma(0.0) // dont do any random oddness
	{
		m_Rand = FastRandom::getInstance(); // Initialise random number generator

		verifyAngleDefs();

		// Conformer definitions
		m_CurrentIndex = 0; // we begin at the 1st conformer, indexed from 0
		m_CurrentConformer.resize(segLength(),conformer_type_fail);

		// Store our internal 'idealised' segment (idealised bonds and angles, a chain break will be produced)
		SegmentDef &internalDef = getSegDef(); // The section we are working on
		internalDef.storeCache(); // cache the current positions in the temporary cache container
		m_Idealised.setTarget(internalDef.getWorkSpace(),internalDef); // Set the idealised position store to point at the segment
		performIdealisation(); // Idealise the cordinates in the workspace - this will result in a single chain break at the end of the segment
		m_Idealised.store(); // Store this idealised state in its special cache for use later as required
		internalDef.revertCache(); // revert the original positions from the temporary cache

		// Report the native confomer to the tra file .. i.e. the one we "want" to find ... if known.
		m_ClosestToNative.setTarget(internalDef.getWorkSpace(),internalDef);
		if( internalDef.isNativeKnown() && !assessNearestToNativeConformer() )
		{
			THROW(ProcedureException,"The conformer builder cannot perform native analysis, the segment does not have a defied native state!");
		}
	}

	void ConfBuilderBase::verifyAngleDefs()
	{
		// Validate that we have the required residues in the angleset library
		const Sequence::BioSequence& seq = getSegDef().getSequence();
		for( size_t i = 0; i < seq.size(); i++ )
		{
			if( !m_AngleSet->getIsDefined( seq[i] ) )
			{
				std::stringstream buff;
				buff << "The residue at position '" << i << "' in the given sequence is not deined in the angleset: '" << seq[i] << "'" << endl;
				THROW(ProcedureException,buff.str());
			}
		}
	}

	void ConfBuilderBase::performIdealisation()
	{
 		printf("Performing idealisation procedure on the conformer...\n");

		size_t index = 0;
		std::vector< double > angles;
		angles.resize(m_RotatePhi_AA.size()+m_RotatePsi_AA.size()+m_RotateOmega_AA.size(),DBL_MAX);
		for( size_t i = 0; i < m_RotatePhi_AA.size(); i++ )
		{
			if( m_RotatePhi_AA[i].isValid() )
				angles[index] = m_RotatePhi_AA[i].getCurrentTorsionAngle();
			index++;
		}
		for( size_t i = 0; i < m_RotatePsi_AA.size(); i++ )
		{
			if( m_RotatePsi_AA[i].isValid() )
				angles[index] = m_RotatePsi_AA[i].getCurrentTorsionAngle();
			index++;
		}
		for( size_t i = 0; i < m_RotateOmega_AA.size(); i++ )
		{
			if( m_RotateOmega_AA[i].isValid() )
				angles[index] = m_RotateOmega_AA[i].getCurrentTorsionAngle();
			index++;
		}

		size_t startRes = getSegDef().getStartResIndex();
		size_t endRes = getSegDef().getEndResIndex();

		SegBreakType breakType = getSegDef().getBreakType();
		if( breakType == SegBreakEnd )
		{
			// Break at end
			// ||--------- // ||
			if( !rebuildExtdFromStart(getWSpace(),startRes,endRes) )
			{
				THROW(ProcedureException,"Idealisation Fail");
			}
		}
		else if( breakType == SegBreakStart )
		{
			// Break at start
			// || // ---------||
			if( !rebuildExtdFromEnd(getWSpace(),startRes,endRes) )
			{
				THROW(ProcedureException,"Idealisation Fail");
			}
		}
		else if( breakType == SegBreakCentre )
		{
			// Break in middle
			// || ----- // ----||
			size_t centreRes = getSegDef().getCentreResIndex();
			if( !rebuildExtdFromStart(getWSpace(),startRes,centreRes) )
			{
				THROW(ProcedureException,"Idealisation Fail");
			}
			if( !rebuildExtdFromEnd(getWSpace(),centreRes+1,endRes) )
			{
				THROW(ProcedureException,"Idealisation Fail");
			}
		}
		else
		{
			THROW(NotImplementedException,"Unknown SegBreakType encountered!");
		}

		index = 0;
		for( size_t i = 0; i < m_RotatePhi_AA.size(); i++ )
		{
			double d = angles[index++];
			if( d != DBL_MAX )
				m_RotatePhi_AA[i].performRotation( d );
		}
		for( size_t i = 0; i < m_RotatePsi_AA.size(); i++ )
		{
			double d = angles[index++];
			if( d != DBL_MAX )
				m_RotatePsi_AA[i].performRotation( d );
		}
		for( size_t i = 0; i < m_RotateOmega_AA.size(); i++ )
		{
			double d = angles[index++];
			if( d != DBL_MAX )
				m_RotateOmega_AA[i].performRotation( d );
		}

		return;
	}

	void ConfBuilderBase::MinimiseIdealised()
	{
		if( !m_IdealisedMinimised )
		{
			// 1) Store current
			SegmentDef &internalDef = getSegDef(); // The section we are working on
			internalDef.storeCache(); // cache the current positions in the temporary cache container

			// 2) Revert standard ideal structure
			m_Idealised.revert();

			// 3) Launch the minimisation
			if( OutputLevel )
			{
				printf("Performing bond angle and torsion minimisation on the idealised conformer\nfor true forcefield compliance...\n(the default ff-file coordinates may be mildly out)\n");
			}

			// Proxies
			size_t startRes = getSegDef().getStartResIndex();
			size_t endRes = getSegDef().getEndResIndex();
			size_t centreRes = getSegDef().getCentreResIndex();

			SegBreakType segBreakType = getSegDef().getBreakType();

			int startBreak = INT_MAX;
			int endBreak = INT_MAX;

			if( segBreakType == SegBreakStart )
			{
				startBreak = getWSpace().res[startRes-1].iC;
				endBreak = getWSpace().res[startRes].iN;
			}
			else if( segBreakType == SegBreakEnd )
			{
				startBreak = getWSpace().res[endRes].iC;
				endBreak = getWSpace().res[endRes+1].iN;
			}
			else if( segBreakType == SegBreakCentre )
			{
				startBreak = getWSpace().res[centreRes].iC;
				endBreak = getWSpace().res[centreRes+1].iN;
			}
			else
			{
				THROW( CodeException, "Unknown SegBreakType encountered!");
			}

			if( OutputLevel ) 
				printf("Calling forcefield setup.\n");

			// Forcefield config
			Physics::Forcefield ff1 = Physics::Forcefield(getWSpace());
			Physics::FF_BreakableBonded* bonds = new Physics::FF_BreakableBonded(getWSpace());
			bonds->OutputLevel = (Verbosity::Type)((int)OutputLevel - 1);
			bonds->createBreak( startBreak, endBreak );
			ff1.addWithOwnership( bonds ) ;

			if( OutputLevel ) 
				printf("Invoking bonded component minimisation.\n");

			// Protocol config
			Protocol::Minimisation minim( ff1, PickResidueRange(getWSpace(), startRes, endRes ) );
			minim.OutputLevel = Verbosity::Silent;
			minim.Steps = 1500; // *maximim*, but will be less than this as 'SlopeCutoff' is defined
			minim.SlopeCutoff = 0.01 / (Physics::PhysicsConst::J2kcal * Physics::PhysicsConst::Na);
			minim.UpdateNList = -1; // Disable; whats the point? We are using bonded-only forces - see ff1 above!
			//minim.UpdateScr = 10;
			//minim.UpdateTra = 10;
			minim.run();
			if( OutputLevel ) 
				printf("Done!\n");

			// We have done it!
			m_IdealisedMinimised = true;

			// 4) Store as the new ideal
			m_Idealised.store(); // Store this idealised state in its special cache for use later as required

			// 5) Revert to the starting structure
			internalDef.revertCache(); // revert the original positions from the temporary cache
		}
	}

	bool ConfBuilderBase::assessNearestToNativeConformer()
	{
		// fill m_CurrentConformer with the "native" conformer
		// obtain this lowest aRMS value
		// print this to the screen and append to the BristolTrajectoryFormat file ...
		// also retrieve the native bin info and fill m_NativeNumChar and m_NativeBinChar

		getSegDef().storeCache(); // backup current structure
		m_Idealised.revert(); // just for sanity in-case something else has tampered with it beforehand ...

		const Sequence::BioSequence& seq = getSegDef().getSequence();

		double origPhi, origPsi, origOmega;
		std::vector<double> closestPhi( segLength() );   // we need to cache these desired values as we want to know the
		std::vector<double> closestPsi( segLength() );   // psi of the final loop residue in the conformer. Thats because it
		std::vector<double> closestOmega( segLength() ); // becomes corupted following rotation...
		double difference;

		m_NativeConformer.resize( segLength() );
		m_NativeNumChar.resize( segLength() );
		m_NativeBinChar.resize( segLength() );

		m_BestARms = 0.0; // init to 0.0
		for( size_t i = 0; i < segLength(); i++ )
		{
			origPhi = m_RotatePhi_BB[i].getCurrentTorsionAngle();
			origPsi = m_RotatePsi_BB[i].getCurrentTorsionAngle();
			origOmega = m_RotateOmega_BB[i].getCurrentTorsionAngle();

			int nativeIndex = m_AngleSet->getClosestAngleGroup(
				seq[i],
				origPhi, closestPhi[i],
				origPsi, closestPsi[i],
				origOmega, closestOmega[i]
			);

			m_NativeConformer[i] = nativeIndex;
			m_NativeNumChar[i] = getScreenChar( nativeIndex );
			m_NativeBinChar[i] = m_AngleSet->getBackboneTorsionSet(seq[i]).getAngleClass(nativeIndex);

			// append the squared deviations to the arms local variable
			difference = closestPhi[i] - origPhi;
			if( difference < 0.0 ) difference = -difference;
			if( difference > Maths::MathConst::PI ) difference = Maths::MathConst::TwoPI - difference;
			m_BestARms += sqr(difference);

			difference = closestPsi[i] - origPsi;
			if( difference < 0.0 ) difference = -difference;
			if( difference > Maths::MathConst::PI ) difference = Maths::MathConst::TwoPI - difference;
			m_BestARms += sqr(difference);

			difference = closestOmega[i] - origOmega;
			if( difference < 0.0 ) difference = -difference;
			if( difference > Maths::MathConst::PI ) difference = Maths::MathConst::TwoPI - difference;
			m_BestARms += sqr(difference);
		}
		m_BestARms = sqrt( m_BestARms / (3.0* (double)segLength()));

		// now apply the rotations ...
		for( size_t i = 0; i < segLength(); i++ )
		{
			// all atom
			//m_RotatePhi_AA  [i].performRotation( closestPhi[i]   );
			//m_RotatePsi_AA  [i].performRotation( closestPsi[i]   );
			//m_RotateOmega_AA[i].performRotation( closestOmega[i] );

			// backbone only
			SegCoreBBOnlyTC::performRotation( m_RotatePhi_BB[i], closestPhi[i] );
			SegCoreBBOnlyTC::performRotation( m_RotatePsi_BB[i], closestPsi[i] );
			SegCoreBBOnlyTC::performRotation( m_RotateOmega_BB[i], closestOmega[i] );
		}
		// ... and so its now transformed to the "native" conformer.

		// do some user reporting.
		if( getSegDef().hasTraComments() )
		{
			IO::BTF_Block_Comment& cmnt = *getSegDef().getTraComments();
			cmnt.setFormat("The closest conformer '%s'(Bin:%s) in terms of aRMS, with an aRMS (Phi,Psi,Omg) of %lf")
				(m_NativeNumChar)
				(m_NativeBinChar)
				(m_BestARms);
		}

		m_ClosestToNative.store(); // store this 'native' state in case we want it later

		getSegDef().revertCache(); // restore backed up conformer

		return true;
	}

	bool ConfBuilderBase::isCurrentConformerNative()
	{
		// The native conformer refers to the original library angleset.
		// the 'current' conformer refers to the current angle sub-set dependent on the current conformer selection mode.
		// Hence the lookup call to getAngleSetID().
		if( !getSegDef().isNativeKnown() ) return false;
		for( size_t i = 0; i < segLength(); i++ )
		{
			int conf = m_BackboneTorsionSets[i].getAngleSetID(m_CurrentConformer[i]);
			if( conf != m_NativeConformer[i] )
			{
				return false;
			}
		}
		return true;
	}

	void ConfBuilderBase::addEntireAnglesetLibrary()
	{
		// Most ConformerBuilder modes require the addition of all possible library angles.
		m_BackboneTorsionSets.clear();

		const Sequence::BioSequence seq = getSegDef().getSequence();
		for( size_t i = 0; i < segLength(); i++ )
		{
			// Obtain the correct library set.
			const BackboneTorsionLibrary& resAngle = m_AngleSet->getBackboneTorsionSet( seq[i] );

			BackboneTorsionSubSet set;
			set.init(resAngle);
			for( size_t j = 0; j < resAngle.size(); j++ )
			{
				// Add* All* the angle subsets for random perturbation
				set.addAngleGroup( j );
			}
			set.finalise();
			m_BackboneTorsionSets.push_back( set );
		}
	}

	void ConfBuilderBase::cloneBackboneTorsionSubSet( ConfBuilderBase& _Builder  )
	{
		m_BackboneTorsionSets = _Builder.m_BackboneTorsionSets;
	}

	std::string ConfBuilderBase::sprintCurrentConformer()
	{
		// The internal conformer is indexed for the number of angles in the subset given to each residue, however
		// when reporting to the user, we want to give the global conformer present in the angleset file....
		// print the global conformer as defined by tha angles in the angleset file...
		std::string conf;
		conf.resize(segLength());
		for( size_t i = 0; i < segLength(); i++ )
		{
			conf[i] = getScreenChar( m_BackboneTorsionSets[i].getAngleSetID( m_CurrentConformer[i] ) );
		}
		return conf;
	}

	void ConfBuilderBase::printCurrentConformer()
	{
		std::cout << sprintCurrentConformer();
	}

	std::string ConfBuilderBase::sprintCurrentInternalConformer()
	{
		// print the internal conformer of the current angle-sub-set ...
		std::string conf;
		conf.resize(segLength());
		for( size_t i = 0; i < segLength(); i++ )
		{
			conf[i] = getScreenChar( m_CurrentConformer[i] );
		}
		return conf;
	}

	void ConfBuilderBase::printCurrentInternalConformer()
	{
		// print the internal conformer of the current angle-sub-set ...
		std::cout << sprintCurrentInternalConformer() << " :INTERNAL CONFORMER";
	}

	void ConfBuilderBase::printAngleSubset()const
	{
		printf("ConformerBuilder AngleSubSet Dataprint:\n");
		for( size_t i = 0; i < segLength(); i++ )
		{
			m_BackboneTorsionSets[i].info();
		}
	}

	void ConfBuilderBase::info()const
	{
		const Sequence::BioSequence seq = getSegDef().getSequence();
		std::string singleLetterSeq = seq.printToStringSingle();

		printf("Conformer Info:\n");

		if( getSegDef().isNativeKnown() )
		{
			printf("\tWhole loop conformation is known!\n");
		}
		else if( getSegDef().isNativeBackboneKnown() )
		{
			printf("\tLoop backbone-only conformation is known!\n");
		}
		else
		{
			printf("\tNone of the loop conformation is known!\n");
		}

		printf( "\t%s: Seg Sequence\n\t", singleLetterSeq.c_str() );
		for( size_t i = 0; i < segLength(); i++ )
		{
			printf("%d", m_BackboneTorsionSets[i].libSize() );
		}
		printf(": Total number of angles in global anglset library\n\t");
		for( size_t i = 0; i < segLength(); i++ )
		{
			printf("%d", m_BackboneTorsionSets[i].size() );
		}
		printf(": Number of angles in use in the conformer subset\n");
		if( getSegDef().isNativeKnown() )
		{
			printf("\t%s: Native Query Num\n", m_NativeNumChar.c_str());
			printf("\t%s: Native Query Bin\n", m_NativeBinChar.c_str());
			printf("\t%lf: Native aRMS (Phi,Psi,Omg)\n", m_BestARms);
		}
	}

	void ConfBuilderBase::setBuildMode( ConformerBuildMode mode )
	{
		if( mode == m_BuildMode ) return;
		m_BuildMode = mode; ///< assign the new build mode
		m_Idealised.revert(); ///< put all the atoms back in the right places
		applyWholeConformer(); ///< re-apply the current conformer
	}

	void ConfBuilderBase::applyWholeConformer()
	{
		for( size_t i = 0; i < segLength(); i++ )
		{
			performRotation( i );
		}
	}

	void ConfBuilderBase::performRotation( int conformerIndex )
	{
		conformer_type conformerID = m_CurrentConformer[conformerIndex];

		double noisePhi = 0.0;
		double noisePsi = 0.0;
		if( NoiseSigma > 0.0 )
		{
			noisePhi = m_Rand->nextNormal( NoiseSigma );
			noisePsi = m_Rand->nextNormal( NoiseSigma );
		}

		if( m_BuildMode == Backbone )
		{
			// code ONLY to move the main backbone atoms leaving sidechains and BB hydrogens in place.
			// this code is intended soley for high speed during initial crude conformer evaluation.

			SegCoreBBOnlyTC::performRotation(
				m_RotatePhi_BB[conformerIndex],
				m_BackboneTorsionSets[conformerIndex].getPhi( conformerID ) + noisePhi
				); // phi

			SegCoreBBOnlyTC::performRotation(
				m_RotatePsi_BB[conformerIndex],
				m_BackboneTorsionSets[conformerIndex].getPsi( conformerID ) + noisePsi
				); // psi

			SegCoreBBOnlyTC::performRotation(
				m_RotateOmega_BB[conformerIndex],
				m_BackboneTorsionSets[conformerIndex].getOmega( conformerID )
				); // omega
		}
		else // we are in all atom mode
		{
			// This code changes the torsions and moves all atoms ....
			RotationDefinition_AllAtom* rotDef;

			// phi
			rotDef = &m_RotatePhi_AA[conformerIndex];
			if( rotDef->isValid() )
			{
				rotDef->performRotation( m_BackboneTorsionSets[conformerIndex].getPhi( conformerID ) + noisePhi );
			}

			// psi
			rotDef = &m_RotatePsi_AA[conformerIndex];
			if( rotDef->isValid() )
			{
				rotDef->performRotation( m_BackboneTorsionSets[conformerIndex].getPsi( conformerID ) + noisePsi );
			}

			// omega
			rotDef = &m_RotateOmega_AA[conformerIndex];
			if( rotDef->isValid() )
			{
				rotDef->performRotation( m_BackboneTorsionSets[conformerIndex].getOmega( conformerID ) );
			}
		}
	}

	// --------------------------------------------------------
	// ConfBuilderBase_Descriptor
	// --------------------------------------------------------

	ConfBuilderBase_Descriptor::ConfBuilderBase_Descriptor(
		const Library::AngleSet& angleSet,
		SegmentDef& segdef
		):
		// Initialise base classes
		ConfBuilderBase(angleSet, segdef),
		SegmentDefUser(segdef),
		// Initialise members
		m_PossibleConformers(0)
	{
		setDescriptor("*");
	}

	void ConfBuilderBase_Descriptor::infoDescriptor() const
	{
		printf("\t%s: Conformational descriptor\n",m_Descriptor.c_str());
		printf("\t%lu: Conformer permutations in total\n",m_PossibleConformers);
	}

	void ConfBuilderBase_Descriptor::setDescriptor( DescriptorType _mode )
	{
		// Reallocation

		switch( _mode )
		{
		case AllResidues:
			{
				m_Descriptor.resize(segLength(),'?');
			}
		case NativeBin:
			{
				m_Descriptor.resize(segLength(),'N');
			}
		case NativeClosest:
			{
				m_Descriptor.resize(segLength(),'C');
			}
		default:
			THROW(CodeException,"Unknown type");
		}

		initDescriptor();
	}

	void ConfBuilderBase_Descriptor::setDescriptor( const std::string& _confDescriptor )
	{
		m_Descriptor = _confDescriptor; // make a local manipulatable copy

		if( m_Descriptor.size() == 0 )
		{
			THROW(ArgumentException,"The length of the conformer descriptor is zero!");
		}
		else if( m_Descriptor.size() == 1 && m_Descriptor[0] == '*' )
		{
			m_Descriptor.clear();
			m_Descriptor.resize(segLength(),'?');
			// This is now a conformer the same length as the internal segment definition
		}
		else if( m_Descriptor.size() != segLength() )
		{
			THROW(ArgumentException,"The length of the conformer descriptor and the segment length do not match!");
		}

		initDescriptor();
	}

	void ConfBuilderBase_Descriptor::initDescriptor()
	{
		validateDescriptor(); // check it makes sense in the current context.
		obtainBBSet(); // now that we have confirmed our descriptor, obtain the angles associaed with that descriptor.
		calcTotalAngles(); // How many angles does this give us?
	}

	void ConfBuilderBase_Descriptor::calcTotalAngles()
	{
		int error = 0;
		m_PossibleConformers = m_BackboneTorsionSets[0].size();
		if( m_PossibleConformers == 0 ) error++;
		for( size_t i = 1; i < segLength(); i++ )
		{
			size_t size = m_BackboneTorsionSets[i].size();
			if( size == 0 ) error++;
			m_PossibleConformers *= size;
		}
		if( error != 0 )
		{
			StringBuilder sb;
			sb.setFormat("There must be at least one angle possibility per resiude. %d of %d angles have 0.")(error)(m_BackboneTorsionSets.size());
			THROW( ProcedureException, sb.toString() );
		}
	}

	void ConfBuilderBase_Descriptor::validateDescriptor()
	{
		bool usingNative = false;
		for( size_t i = 0; i < segLength(); i++ )
		{
			m_Descriptor[i] = toupper( m_Descriptor[i] ); // Force upper-case (we want global case insensitivity in input files)
			if( m_Descriptor[i] == 'N' || m_Descriptor[i] == 'C' )
			{
				usingNative = true;
				break;
			}
		}
		if( usingNative && !getSegDef().isNativeKnown() )
		{
			THROW(ArgumentException,"The conformer descriptor requires 'Native' ('N'/'C') knowledge. The native is not defined in the loop definition!");
		}
	}

	void ConfBuilderBase_Descriptor::obtainBBSet()
	{
		// Clean up any existing memory
		m_BackboneTorsionSets.clear();

		// now get the angle subset that is defined by this descriptor
		// Definitions are as follows:
		// ? = any
		// N = "same as native bin" (i.e. any one of of A, B, L, -) as defined by the entry cached in the initial state of the "loopdefinition" class instance
		// Alphanumeric or '-'... Alphanumeric characters define conformer classes. '-' defines all the deviant angles or angles with no defined class.
		// Digits = defined a single angle via the respective position in the angleset residue library

		std::string seq = getSegDef().getSequence().printToStringSingle();

		ASSERT( seq.size() == segLength(), CodeException, "Internal sequence structure mismatch");

		for( size_t i = 0; i < segLength(); i++ )
		{
			const BackboneTorsionLibrary& resAngle = m_AngleSet->getBackboneTorsionSet( seq[i] );

			BackboneTorsionSubSet set;
			set.init( resAngle );

			// Add a subset of the the angle group dependent on the descriptor
			if( '?' == m_Descriptor[i] )
			{
				// Add them all
				for( size_t j = 0; j < resAngle.size(); j++ )
				{
					set.addAngleGroup( j );
				}
			}
			else if ( 'C' == m_Descriptor[i] )
			{
				// We have already called validateDescriptor() - we defo have native knowledge
				if( !set.addAngleGroup( m_NativeConformer[i] ) )
				{
					THROW( ProcedureException, "Error adding angleGroup during 'ConfBuild_FromDescriptor' initialisation!" );
				}
			}
			else if( 'N' == m_Descriptor[i] )
			{
				// We have already called validateDescriptor() - we defo have native knowledge
				// Add the angles in the same angleset bin as the native defined by the loop definition.
				for( size_t j = 0; j < resAngle.size(); j++ )
				{
					if( resAngle.getAngleClass( j ) == m_NativeBinChar[i] )
					{
						if( !set.addAngleGroup( j ) )
						{
							THROW( ProcedureException, "Error adding angleGroup during 'ConfBuild_FromDescriptor' initialisation!" );
						}
					}
				}
			}
			else if( isalpha( m_Descriptor[i] ) || m_Descriptor[i] == '-' ) // valid class descriptors
			{
				for( size_t j = 0; j < resAngle.size(); j++ )
				{
					if( resAngle.getAngleClass( j ) == m_Descriptor[i] )
					{
						if( !set.addAngleGroup( j ) )
						{
							THROW( ProcedureException, "Error adding angleGroup during 'ConfBuild_FromDescriptor' initialisation!" );
						}
					}
				}
				if( set.size() == 0 )
				{
					std::stringstream errorMessage;
					errorMessage << "Could not find an example of descriptor class '"
						<< m_Descriptor[i]
						<< "' in the angleset residue library for position "
						<< i << " in the conformer.";
					THROW(ArgumentException, errorMessage.str());
				}
			}
			else if( isdigit( m_Descriptor[i] ) )
			{
				// NOTE :
				// there will be just one angle in the REDUCED angle set that is created here
				// therefore
				// the internal conformer descriptor m_CurrentConformer[i] is Always == 0
				int obtainIndex = ((int)m_Descriptor[i]) - 48; ///< 48 to convert the ASCII character to an index in the true angleset ...
				if( obtainIndex < 0 || obtainIndex >= (int)set.libSize() )
				{
					std::stringstream errorMessage;
					errorMessage << "The defined descriptor digit of the conformer position " << i << " is not in the range defined by the angleset library.";
					THROW(ArgumentException, errorMessage.str());
				}
				if( !set.addAngleGroup( obtainIndex ) )
				{
					THROW( ProcedureException, "Error adding angleGroup during 'ConfBuild_FromDescriptor' initialisation!" );
				}
			}
			else
			{
				std::stringstream errorMessage;
				errorMessage << "Conformer input index '" << i << "' uses an invalid descriptor character: '" << m_Descriptor[i] << "'";
				THROW(ArgumentException, errorMessage.str());
			}
			set.finalise();

			if( set.size() == 0 )
			{
				std::stringstream errorMessage;
				errorMessage << "Conformer input index '" << i << "'s descriptor character '" << m_Descriptor[i] << "' results in 0 valid angles!";
				THROW(ArgumentException, errorMessage.str());
			}

			// All cool :-D
			m_BackboneTorsionSets.push_back(set);
		}
	}

	// --------------------------------------------------------
	// ConfBuilderBase_Random
	// --------------------------------------------------------
	ConfBuilderBase_Random::ConfBuilderBase_Random( const AngleSet &angleSet, SegmentDef &segdef ) :
		ConfBuilderBase_Descriptor(angleSet, segdef),
		SegmentDefUser(segdef),
		m_RandConformers(0), // default to none
		PropensityWeighting(false)
	{
		addEntireAnglesetLibrary(); // For random selection, we need the entire angleset.
		RandomiseWholeConformer(); // initialise to a random conformer
		reset(); // set index to 0
	}

	void ConfBuilderBase_Random::setRandCount( conformer_count_type randomCount, bool warnOverMax )
	{
		m_RandConformers = randomCount;
		if( warnOverMax && getPermutations() < m_RandConformers )
		{
			std::cout << "WARNING: Requested random count exceeds the maximum possible number of conformer permutations!" << std::endl;
		}
	}

	void ConfBuilderBase_Random::reset()
	{
		m_CurrentIndex = 0; // this will allow the next() call again.
	}

	void ConfBuilderBase_Random::RandomiseWholeConformer()
	{
		for( size_t i = 0; i < segLength(); i++ )
		{
			RandomiseSingleConformerPos( i );
		}
	}

	void ConfBuilderBase_Random::RandomiseSingleConformerPos()
	{
		RandomiseSingleConformerPos( m_Rand->next(0,segLength()) );
	}

	void ConfBuilderBase_Random::RandomiseSingleConformerPos( size_t _changeIndex )
	{
		D_ASSERT( _changeIndex < segLength(), CodeException, "Conformer change index is ouside conformer range!");
		m_ChangeIndex = _changeIndex;

		size_t nConf = m_BackboneTorsionSets[m_ChangeIndex].size();
		if( nConf == 1 ) 
		{
			// Special case, otherwise we will get an infinite loop below!
			m_CurrentConformer[m_ChangeIndex] = 0;
			performRotation( m_ChangeIndex );
			return;
		}

		conformer_type original = m_CurrentConformer[m_ChangeIndex];
		conformer_type changeTo;
		if( PropensityWeighting )
		{
			// We *KNOW*, by definition of the angleet that the propensities sum to 1.0
			while( true )
			{
				double fract = m_Rand->nextDouble();
				double fractSum = 0.0;
				Library::BackboneTorsionSubSet& ang = m_BackboneTorsionSets[m_ChangeIndex];
				changeTo = UCHAR_MAX;
				size_t toThis = ang.size() - 1;
				for( size_t i = 0; i < toThis; i++ )
				{
					fractSum += ang.getPropensity(i);
					if( fractSum > fract )
					{
						changeTo = i;
						break;
					}
				}
				if( changeTo == UCHAR_MAX )
					changeTo = ang.size()-1;
				if( changeTo != original )
					break;
			}
		}
		else
		{
			while( true )
			{
				changeTo = m_Rand->next( 0, nConf );
				if( changeTo != original )
					break;
			}
		}

		m_CurrentConformer[m_ChangeIndex] = changeTo;
		performRotation( m_ChangeIndex );
	}

	// --------------------------------------------------------
	// ConfBuilderBase_Enum
	// --------------------------------------------------------

	ConfBuilderBase_Enum::ConfBuilderBase_Enum( Library::AngleSet& angleSet, SegmentDef& segdef ):
		// Initialise base classes
		ConfBuilderBase_Descriptor(angleSet, segdef),
		SegmentDefUser(segdef),
		// Set the enumeration mode to default to the fastest method
		m_EnumerationMode(NumericAscent)
	{
	}

	void ConfBuilderBase_Enum::ClearSequentialModeMemory()
	{
		// Memory from SequentialAscent Mode
		m_State_MaxIndx.clear();
		m_State_RepeatCounter.clear();
		m_State_Root.clear();
		m_State_Current.clear();
		m_SortTemplate_Current.clear();
		m_SortTemplate_Root.clear();
		m_SubStateTaken.clear();
		m_SubState_HadLast.clear();
	}

	void ConfBuilderBase_Enum::changeEnumerationMode( ConformerEnumMode _Mode )
	{
		reset();
		m_EnumerationMode = _Mode;

		if( _Mode == SequentialAscent )
		{
			m_State_MaxIndx.resize(segLength());
			m_State_Root.resize(segLength());
			m_State_Current.resize(segLength());
			m_SortTemplate_Current.resize(segLength());
			m_SortTemplate_Root.resize(segLength());
			m_SubState_HadLast.resize(segLength());
			m_SubStateTaken.resize(segLength());

			initStates();
		}
		else
		{
			// other modes dont require additional memory.
			ClearSequentialModeMemory();
		}

		return;
	}

	void ConfBuilderBase_Enum::reset()
	{
		for( size_t i = 0; i < segLength(); i++ )
		{
			m_CurrentConformer[i] = 0; ///< initialise as the first angle in the subset defined by the descriptor, it will define a minimum of 1 angle for each position in the conformer
		}
		// both these are used in conformer descriptor incrementation
		m_CurrentIndex = 0; ///< set this
		applyWholeConformer(); ///< re-apply the initial conformer
	}

	bool ConfBuilderBase_Enum::next()
	{
		m_CurrentIndex++; ///< increment

		switch( m_EnumerationMode )
		{
		case NumericAscent:
			{
				if( m_CurrentIndex >= getPermutations() )
				{
					return false; ///< we have reached the end
				}
				else
				{
					incrementDescriptorAndapply(0);
					return true;
				}
			}
		case SequentialAscent:
			{
				// Make each STATE e.g. m_State_Root: 2110
				// And enumerate it, yielding:
				// 2110, 2101, 2011, 1210, 1201, 1120, 1102, 0211, 0121, 0112
				if (nextSubState())
				{
					applySequentialAscentConformer();
					return true; ///< the next valid enumeration from the current state
				}
				else if (nextState())
				{
					initSubState(); ///< the first valid enumberation from the current state
					applySequentialAscentConformer();
					return true;
				}
				else
				{
					return false; ///< no more states
				}
			}
		default:
			{
				THROW(CodeException,"Unknown 'ConformerEnumMode' encountered!");
			}
		}
	}

	void ConfBuilderBase_Enum::incrementDescriptorAndapply( int conformerIndex )
	{
		m_CurrentConformer[ conformerIndex ]++;// increment the current one ...
		// if we have passed the limit, reset and increment the next index ...
		if( m_CurrentConformer[ conformerIndex ] < (conformer_type) m_BackboneTorsionSets[ conformerIndex ].size() )
		{
			performRotation( conformerIndex );
		}
		else
		{
			m_CurrentConformer[ conformerIndex ] = 0; ///< reset the current index and pass-the-buck
			performRotation( conformerIndex );
			incrementDescriptorAndapply( conformerIndex + 1 );
		}
		return;
	}

	void ConfBuilderBase_Enum::applySequentialAscentConformer()
	{
		int changed = 0;
		for( size_t i = 0; i < segLength(); i++ )
		{
			if( m_CurrentConformer[i] != m_State_Current[i] )
			{
				m_CurrentConformer[i] = m_State_Current[i];
				performRotation( i );
				changed++;
			}
		}
		//printf("changed %d!\n",changed);
	}

	void ConfBuilderBase_Enum::initStates()
	{
		// MakeStatesMax // Sort it. It 'simplifies' the enumeration.
		for (size_t i = 0; i < segLength(); i++)
		{
			m_State_MaxIndx[i] = m_BackboneTorsionSets[i].size() - 1;
		}

		// sort them
		bool changes = true;
		conformer_type swap = 0;
		while (changes)
		{
			changes = false;
			for (size_t i = 0; i < segLength(); i++)
			{
				for (size_t j = i + 1; j < segLength(); j++)
				{
					if (m_State_MaxIndx[i] < m_State_MaxIndx[j])
					{
						swap = m_State_MaxIndx[i];
						m_State_MaxIndx[i] = m_State_MaxIndx[j];
						m_State_MaxIndx[j] = swap;
						changes = true;
					}
				}
			}
		}

		m_State_MaxMax = 0;
		for (size_t i = 0; i < segLength(); i++)
		{
			m_SubStateTaken[i] = false;
			m_State_Root[i] = 0;
			m_SortTemplate_Root[i] = i; // 0,0,0,0 should be 0,1,2,3
			if (m_State_MaxIndx[i] > m_State_MaxMax) m_State_MaxMax = m_State_MaxIndx[i];
		}

		m_State_RepeatCounter.resize(m_State_MaxMax + 1); // +1 as MaxMax is an index, not a length
		m_StateResultantIndex = 0;
		for (size_t i = 0; i < m_State_MaxMax; i++)
		{
			m_State_RepeatCounter[i] = 0;
		}

		initSubState();
	}

	bool ConfBuilderBase_Enum::nextState()
	{
		if (m_StateResultantIndex == 0 && (m_State_Root[m_StateResultantIndex] == (m_State_MaxIndx[m_StateResultantIndex])))
		{
			return false; ///< the only inner-loop break condition ...
		}

		bool changed = false;
		if (m_State_Root[m_StateResultantIndex] < m_State_MaxIndx[m_StateResultantIndex])
		{
			m_State_Root[m_StateResultantIndex]++;
			for (size_t i = (size_t)m_StateResultantIndex + 1; i < segLength(); i++)
			{
				m_State_Root[i] = 0;
			}

			//Make SortTemplateRoot:
			for (conformer_type i = 0; i <= m_State_MaxMax; i++)
			{
				m_State_RepeatCounter[i] = 0;
			}
			for (size_t i = 0; i < segLength(); i++)
			{
				m_SortTemplate_Root[i] = m_State_RepeatCounter[m_State_Root[i]]++;
			}

			//printStateDefinition();
			changed = true;
		}

		// find the next index ...
		m_StateResultantIndex++;
		if (m_StateResultantIndex >= (int)segLength()) // should never be larger
		{
			for (m_StateResultantIndex = segLength() - 1; m_StateResultantIndex > 0; m_StateResultantIndex--)
			{
				if (
					(m_State_Root[m_StateResultantIndex] < m_State_Root[m_StateResultantIndex - 1]) &&
					(m_State_Root[m_StateResultantIndex] < m_State_MaxIndx[m_StateResultantIndex])
					) break;
			}
		}

		if (changed)
		{
			return true;
		}
		else
		{
			return nextState();
		}
	}

	void ConfBuilderBase_Enum::initSubState()// with the next valid entry
	{
		// setup the first state from the 'm_State_Root'
		m_SubStatePosition = segLength();
		for (size_t i = 0; i < segLength(); i++)
		{
			m_SubState_HadLast[i] = i;
			m_SubStateTaken[i] = true;
			m_State_Current[i] = m_State_Root[i];
			m_SortTemplate_Current[i] = m_SortTemplate_Root[i];
		}

		// make sure it is valid, they occationally arent, but at least one obviously must be.
		if (!isValid()) nextSubState();
	}

	bool ConfBuilderBase_Enum::nextSubState()
	{
		do
		{
			m_SubStatePosition -= 1;
			while (m_SubStatePosition != segLength())
			{
				if (m_SubState_HadLast[m_SubStatePosition] != -1)
				{
					m_SubStateTaken[m_SubState_HadLast[m_SubStatePosition]] = false;
				}

				bool setCheck = false;
				for (size_t i = m_SubState_HadLast[m_SubStatePosition] + 1; i < segLength(); i++)
				{
					if (!m_SubStateTaken[i])
					{
						m_SubStateTaken[i] = true;
						m_SubState_HadLast[m_SubStatePosition] = i;
						m_State_Current[m_SubStatePosition] = m_State_Root[i];
						m_SortTemplate_Current[m_SubStatePosition] = m_SortTemplate_Root[i];
						setCheck = true;
						break;
					}
				}

				if (setCheck)
				{
					m_SubStatePosition++;
				}
				else
				{
					m_SubState_HadLast[m_SubStatePosition] = -1;
					m_SubStatePosition--;
					if (m_SubStatePosition == -1)
						return false; ///< our break condition, there are no more entries!
				}
			}
		}
		while (!isValid());
		return true;
	}

	bool ConfBuilderBase_Enum::isValid()
	{
		// SortTemplateisValid
		int want;
		for (conformer_type j = 0; j <= m_State_MaxMax; j++)
		{
			want = 0;
			for (size_t i = 0; i < segLength(); i++)
			{
				if (m_State_Current[i] == j)
				{
					if (want != m_SortTemplate_Current[i])
					{
						return false;
					}
					want++;
				}
			}
		}

		// ConformerisValid
		for (size_t i = 0; i < segLength(); i++)
		{
			if (m_State_Current[i] >= m_BackboneTorsionSets[i].size())
			{
				return false;
			}
		}

		return true;
	}

	void ConfBuilderBase_Enum::printStateDefinition()
	{
		printf("STATE: ");
		for (size_t i = 0; i < segLength(); i++)
		{
			printf("%d", m_State_Root[i]);
		}
		printf("\n");
	}


	// --------------------------------------------------------
	// ConfBuilderBase_Buffer
	// --------------------------------------------------------

	ConfBuilderBase_Buffer::ConfBuilderBase_Buffer(const  Library::AngleSet& angleSet, SegmentDef &segdef):
		// Initialise base classes
		ConfBuilderBase(angleSet, segdef),
		SegmentDefUser(segdef),
		// Initialise members
		m_InternalCapacity(0) // 0 (and <0) actually mean 'no limit' - it makes no sense to have no capacity - and 0 makes a nice null value.
	{
		setInfinateCapacity();
	}

	void ConfBuilderBase_Buffer::setMaxCapacity( size_t maxBufferCapacity )
	{
		m_InternalCapacity = maxBufferCapacity; // Set the internal value

		// initialise the buffer
		if( m_InternalCapacity != 0 )
		{
			m_Conformers.reserve( maxBufferCapacity ); //reserve the correct amount of memory
		}
		else
		{
			return; // Setting to 0 means that we have infinate capacity, just return
		}

		// Deleted memory
		if( m_Conformers.size() > maxBufferCapacity )
		{
			m_Conformers.erase(m_Conformers.begin()+maxBufferCapacity,m_Conformers.end());
		}

		// Special tool call to remove memory from the internal buffer
		TrimVectorCapacity(m_Conformers);
	}

	bool ConfBuilderBase_Buffer::next()
	{
		m_CurrentIndex++; ///< increment
		if( m_CurrentIndex >= m_Conformers.size() ) return false; ///< we have reached the end

		rotateFromBuffer(m_CurrentIndex);
		return true;
	}

	void ConfBuilderBase_Buffer::reset()
	{
		m_CurrentIndex = 0;
	}

	bool ConfBuilderBase_Buffer::validateConformer( const Conformer& newBuffer )
	{
		for( size_t i = 0; i < segLength(); i++ )
		{
			if( newBuffer[i] >= m_BackboneTorsionSets[i].size() )
			{
				return false;
			}
		}
		return true;
	}

	bool ConfBuilderBase_Buffer::push_back( const Conformer& _Conformer )
	{
		// m_InternalCapacity can be <= 0 to turn off capacity validation
		if( (m_InternalCapacity > 0) && (size() == m_InternalCapacity) )
		{
			printf("WARNING: Conformer buffer overflow. New entry ignored!");
			return false; ///< internal capacity overflow. (We want this in case we want to limit the number of conformers passed to the next stage.
		}
		else
		{
			// check the new buffer for validity before we add it to the list!
			if( !validateConformer( _Conformer ) )
			{
				printf( "WARNING: Call to append() gave an invalid Conformer. New entry ignored!" );
				return false;
			}

			m_Conformers.push_back( _Conformer ); ///< add the buffer to the current list

			return true;
		}
	}

	bool ConfBuilderBase_Buffer::rotateFromBuffer( conformer_count_type conformerIndex )// set the ConformersCurrent state to the conformer in the buffer
	{
		if( conformerIndex >= m_Conformers.size() ) // conformerIndex < 0 cannot be because conformer_count_type is unsigned
		{
			THROW(OutOfRangeException,"Index out of bounds");
		}

		// Unsafe?
		printf("WARNING - check this is working ;-D");
		memcpy(&m_CurrentConformer[0], &m_Conformers[conformerIndex][0], sizeof(conformer_type) * segLength());
		// Safe
		//m_CurrentConformer = m_Conformers[i];

		applyWholeConformer();
		return true;
	}

	void ConfBuilderBase_Buffer::printBuffer()
	{
		std::string seqString = getSegDef().getSequence().printToStringSingle();
		printf("\nConformerBuffer print: (These are the indexes in the ConformerBuilder AngleSubset NOT the original AngleSet.)\n");
		printf("%s :Sequence\n",seqString.c_str());
		for( size_t i = 0; i < m_Conformers.size(); i++ )
		{
			Conformer& conf = m_Conformers[i];
			for( size_t j = 0; j < conf.size(); j++ )
			{
				printf("%c",getScreenChar(conf[j]));
			}
			printf("\n");
		}
	}
}

