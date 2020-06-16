#include "global.h"    

#include "tools/io.h"
#include "library/angleset.h"
#include "workspace/workspace.h" 
#include "workspace/neighbourlist.h"
#include "forcefields/ffbonded.h"
#include "forcefields/breakablebonded.h"
#include "protocols/torsionalminimisation.h"
#include "protocols/minimise.h"

#include "filters/surface.h"
#include "filters/segdistfan.h"
#include "system/segbuilder.h"

#include "protocols/arcus_base.h"  

using namespace Manipulator;

namespace Protocol
{
	ArcusBase::ArcusBase( Physics::Forcefield& _ffs, Physics::Forcefield & _ff, const Library::AngleSet& _as ):
		ProtocolBase(_ff),
		m_AngleSet( &_as ),
		ffs(&_ffs),
		m_SplitBuilder(false),
		SetRefStatePostInitialRebuild(true)
	{
		ExtraReport = false;
		FlagRebuilt = false;
		UtiliseIdealised = true;
		MinimiseIdealised = true;

		m_AutoRegions = true;
		Step = 0;
		Steps = 1;

		addFFBreaks();
	}

	void ArcusBase::info() const
	{
		printf("\n'ArcusBase' Info: \n");
		printf("---------------\n");
		printf("  Generate : %d Models \n", Steps); 
	}

	void ArcusBase::buidAutoDetect()
	{
		m_AutoRegions = true;
		clear();
	}

	void ArcusBase::buildAddExplicit( const PickResidueRange& _Range )
	{
		m_AutoRegions = false;
		m_Regions.push_back( SegmentDef( getWSpace(), _Range.getStartResIndex(), _Range.getEndResIndex(), SegBreakCentre ) );
	}

	void ArcusBase::clear()
	{
		if( m_AutoRegions ) 
		{
			m_Regions.clear();
		}
		m_ConfBuilder.clear();
		m_ConfBuilderSeg.clear();
		m_Filters.clear();
	}
	
	int ArcusBase::run()
	{
		if(Steps < 0)
		{
			throw(ArgumentException("'Steps' determines the number of candidate models that are generated. This must be greater than 0."));
		}

		if(ensureFFSetup()!=0) 
			return -1;
		
		// Initialisation Phase
		int init = initialise();
		if( init <= 0 ) // returns 0 if there are no regions to manipulate
			return init;

		ASSERT( 
			m_ConfBuilder.size() == m_Filters.size(), 
			CodeException, "Internal paired-array size error");

		int steps = runcore();

		// If we have FlagRebuilt flagged, do some setRebuildRequired() 
		callFlagRebuiltIf();

		// Do some housekeeping
		clear();

		return steps;
	}

	bool ArcusBase::verifyNoRegionClash()
	{
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			for( size_t j = i+1; j < m_Regions.size(); j++ )
			{
				size_t rIStart = m_Regions[i].getStartAtomIndex();
				size_t rJStart = m_Regions[j].getStartAtomIndex();
				size_t rIEnd = m_Regions[i].getEndAtomIndex();
				size_t rJEnd = m_Regions[j].getEndAtomIndex();
				if( rJStart <= rIEnd && rJStart >= rIStart ) 
					return false; // rJStart cannot lie within the range of rI
				if( rJEnd <= rIEnd && rJEnd >= rIStart ) 
					return false; // rJEnd cannot lie within the range of rI
			}
		}
		return true;
	}

	void ArcusBase::configLevel1Filters( size_t i )
	{
		// Make the relevent filters for each segment
		m_Filters.push_back( FilterContainer() );

		// 1) Join Filter
		OmegaGroupFilter* joinFilter = new OmegaGroupFilter();
		joinFilter->setTo( getWSpace(), m_Regions[i].getBreakResIndex(), 3.0 );
		m_Filters[i].addWithOwnership( joinFilter );

		// 2) SegmentDistanceFilter: Ensure join is possible
		if( SegFanFilename.size() > 0 )
		{
			bool reverse = i % 2 != 0;
			ASSERT( IO::fileExists( SegFanFilename ), IOException, "Cannot find seg-dat file");
			SegmentDistanceFilter* distFilter = new SegmentDistanceFilter(SegFanFilename);
			distFilter->initialise( m_Regions[i], reverse );
			m_Filters[i].addWithOwnership( distFilter );
		}

		// 3) Surface clash filter
		SurfacePicker surfaceBox;
		surfaceBox.calcBoundingBox( m_Regions[i] );
		ClashGridFilter* filter = new ClashGridFilter(); 
		filter->setForPicker( Pick_AND(PickResidueRange(m_Regions[i]),PickCoreBackbone()) );
		filter->setAgainstPicker( Pick_AND( Pick_NOT(PickResidueRange(m_Regions[i])),surfaceBox) );
		filter->setMolecule( getWSpace() );
		m_Filters[i].addWithOwnership( filter );
	}

	bool ArcusBase::initialBuild()
	{
		// Revert to an uninitialised state
		clear();

		if( m_AutoRegions )
		{
			// Autodetection is enabled
			SegmentRebuilder builder(false);
			builder.DownstreamFullSectionRebuild = true; // dont mark isRebuildRequired false for atoms in build sections
			if( !builder.invokeBuild( getWSpace(), OutputLevel) ) 
				return false; // rebuild missing sections and atoms

			const std::vector<PickResidueRange>& segs = builder.getPreviousBuildSections();
			for( size_t i = 0; i < segs.size(); i++ )
			{
				// Upgrade the PickResidueRange to a full SegDef
				m_Regions.push_back( 
					SegmentDef( 
					getWSpace(), 
					segs[i].getStartResIndex(), 
					segs[i].getEndResIndex(), 
					SegBreakCentre ) );
			}
		}
		else
		{
			ASSERT( m_Regions.size() > 0, CodeException, "Breakdown in code assumption: if 'm_AutoRegions' is flagged there should always be data in 'm_Regions'.");
			// buildExplicit() has therefore been called at some point.

			MissingAtomRebuilder buildMinor;
			if( !buildMinor.invokeBuild( getWSpace(), OutputLevel) ) 
				return false; // rebuild missing sections and atoms

			// Mark sections which have been built as not-built. Only non-builder atoms should be marked
			// as rebuilt by the above builder... (i.e. only non-loop H-atoms elsewhere in the structure)
			setFlagRebuilt(false);

			if( !verifyNoRegionClash() )
			{
				if( OutputLevel ) 
					printf("'ArcusBase' initial build failed!\n");
				return false; // Error condition
			}
		}

		if( SetRefStatePostInitialRebuild )
			getWSpace().setAsReference(); // for crms operations

		m_ConfBuilderSeg.reserve( m_Regions.size() * 2 ); // otherwise buffer could be reallocated, invalidating pointers ahead!!
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			// Create the required array of conformer builders
			if( m_SplitBuilder && m_Regions[i].getBreakType() == SegBreakCentre )
			{
				size_t ii = m_ConfBuilderSeg.size();
				m_ConfBuilderSeg.push_back( SegmentDef( getWSpace(), m_Regions[i].getStartResIndex(), m_Regions[i].getCentreResIndex(), SegBreakEnd ) );
				m_ConfBuilder.push_back( ConfBuild_RandomSingle( *m_AngleSet, m_ConfBuilderSeg[ii] ) );
				
				ii = m_ConfBuilderSeg.size();
				m_ConfBuilderSeg.push_back( SegmentDef( getWSpace(), m_Regions[i].getCentreResIndex()+1, m_Regions[i].getEndResIndex(), SegBreakStart ) );
				m_ConfBuilder.push_back( ConfBuild_RandomSingle( *m_AngleSet, m_ConfBuilderSeg[ii] ) );
			}
			else
			{
				// Make a conformerbuilder for each segment
				m_ConfBuilder.push_back( ConfBuild_RandomSingle( *m_AngleSet, m_Regions[i] ) );
			}
		}

		for( size_t i = 0; i < m_ConfBuilder.size(); i++ )
		{
			m_ConfBuilder[i].setRandCount( m_ConfBuilder[i].getPermutations() );
			if( UtiliseIdealised )
			{
				if( MinimiseIdealised )
				{
					m_ConfBuilder[i].MinimiseIdealised();
				}
				m_ConfBuilder[i].revertIdealised(); // we want to perform the buildup on the idealised structure
			}
			configLevel1Filters(i);
		}

		return true;
	}

	int ArcusBase::addFFBreaks()
	{
		// Get all forcefield components that support the 'Physics::IBondBreakable' *interface* (Abstract Base Class)
		Physics::obtainFFComponents(*ffs,m_BondedFFS);
		return addFFBreaks( m_BondedFFS );

		// Physics::obtainFFComponents(_ff,m_BondedFF); <- most likely this is a baaad plan as default
		//addFFBreaks( m_BondedFF );
	}

	int ArcusBase::addFFBreaks( std::vector<Physics::IBondBreakable*>& _bondedComponents )
	{
		// None of the forcefields require notification of a break.
		if( _bondedComponents.size() == 0 ) 
			return 0;

		int breaks = 0;
		for( size_t j = 0; j < _bondedComponents.size(); j++ )
		{
			_bondedComponents[j]->clearBreaks();
		}

		// Add the relevent break definitions
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			SegBreakType breakType = m_Regions[i].getBreakType();
			if( breakType == SegBreakStart )
			{
				int startIndex = m_Regions[i].getStartResIndex();
				if( !getWSpace().isMolStartResIndex(startIndex) ) 
				{
					// If its a start residue there is no point in adding a break - it IS a terminus already!
					int iC = getWSpace().res[startIndex-1].iC;
					int iN = getWSpace().res[startIndex].iN;
					if( iC < 0 ) return -1; // Cant find atom!
					if( iN < 0 ) return -1; // Cant find atom!
					breaks++;
					for( size_t j = 0; j < _bondedComponents.size(); j++ )
					{
						_bondedComponents[j]->createBreak(iC,iN);
					}
				}
			}
			else if( breakType == SegBreakEnd )
			{
				int endIndex = m_Regions[i].getEndResIndex();
				if( !getWSpace().isMolEndResIndex(endIndex) ) 
				{
					// If its an end residue there is no point in adding a break - it IS a terminus already!
					int iC = getWSpace().res[endIndex].iC;
					int iN = getWSpace().res[endIndex+1].iN;
					if( iC < 0 ) return -1; // Cant find atom!
					if( iN < 0 ) return -1; // Cant find atom!
					breaks++;
					for( size_t j = 0; j <_bondedComponents.size(); j++ )
					{
						_bondedComponents[j]->createBreak(iC,iN);
					}
				}
			}
			else if( breakType == SegBreakCentre )
			{
				int centreIndex = m_Regions[i].getCentreResIndex();
				// If its an end residue there is no point in adding a break - it IS a terminus already!
				int iC = getWSpace().res[centreIndex].iC;
				int iN = getWSpace().res[centreIndex+1].iN;
				if( iC < 0 ) return -1; // Cant find atom!
				if( iN < 0 ) return -1; // Cant find atom!
				breaks++;
				for( size_t j = 0; j < _bondedComponents.size(); j++ )
				{
					_bondedComponents[j]->createBreak(iC,iN);
				}
			}
			else
			{
				THROW( NotImplementedException, "Unknown SegBreakType encountered" );
			}
		}

		return breaks;
	}

	int ArcusBase::initialise()
	{
		// 1) Get segment starting conformations
		if( OutputLevel ) printf("Initialising starting conformations for missing sections.\n");
		if( !initialBuild() )
		{
			if( OutputLevel ) printf("'ArcusBase' initial build failed!\n");
			return -1; // Error condition
		}
		if( m_Regions.size() == 0 )
		{
			if( OutputLevel ) printf("'ArcusBase' exiting early: No sections require rebuild!\n");
			// If this is not the case, check that the System created using the standard rebuilder
			// All atoms in missing regions must have isRebuildRequired() flagged.
			return 0; // no steps have been performed.
		}
		if( ExtraReport ) getWSpace().outtra.append();

		// 2) Setup the forcefield in the context of the bond break-list
		if( addFFBreaks() < 0 )
		{
			if( OutputLevel ) printf("'ArcusBase' could not create the required forcefield breaks!\n");
			return -1;
		}
		if( OutputLevel ) printf("\nCalling ff setup.\n");
		if(ensureFFSetup()!=0) return -1; // Error condition

		return m_Regions.size();
	}

	void ArcusBase::setPickedRegions( const WorkSpace& wspace, PickAtomRanges& picker, SegmentDef& _region ) const
	{
		picker.clear();
		setPickedRegions_Core(wspace,picker,_region);
		picker.assertRanges(wspace);
	}

	void ArcusBase::setPickedRegions_Core( const WorkSpace& wspace, PickAtomRanges& picker, SegmentDef& _region ) const
	{
		SegBreakType bt = _region.getBreakType();
		switch( bt )
		{
		case SegBreakStart:
			{
				picker.addRange(_region,true);
				break;
			}
		case SegBreakEnd:
			{
				picker.addRange(_region,false);
				break;
			}
		case SegBreakCentre:
			{
				size_t breakRes = _region.getBreakResIndex();
				size_t startAtomIndex = _region.getStartAtomIndex();
				size_t centreAtomIndexA = wspace.res[breakRes-1].ilast;
				picker.addRange(startAtomIndex,centreAtomIndexA-startAtomIndex+1,false);
				size_t endAtomIndex = _region.getEndAtomIndex();
				size_t centreAtomIndexB = wspace.res[breakRes].ifirst;
				ASSERT( centreAtomIndexA == centreAtomIndexB-1, CodeException, "Oddness happened" );
				picker.addRange(centreAtomIndexB,endAtomIndex-centreAtomIndexB+1,true);
				break;
			}
		default:
			THROW( CodeException, "Unknown SegBreakType encountered");
		}	
	}

	void ArcusBase::setPickedRegions( const WorkSpace& wspace, PickAtomRanges& picker, std::vector<SegmentDef> _regions ) const
	{
		picker.clear();
		for( size_t j = 0; j < _regions.size(); j++ )
		{
			setPickedRegions_Core(wspace,picker,_regions[j]);		
		}
		picker.assertRanges(wspace);
	}

	void ArcusBase::callFlagRebuiltIf()
	{
		if( FlagRebuilt )
		{
			setFlagRebuilt(true);
		}
	}

	void ArcusBase::setFlagRebuilt( bool setBuilt )
	{
		for( size_t i = 0; i < m_Regions.size(); i++ )
		{
			for( size_t j = m_Regions[i].getStartAtomIndex(); j <= m_Regions[i].getEndAtomIndex(); j++ )
			{
				getWSpace().getAtom(j).setRebuildRequired(!setBuilt);
			}
		}
	}

	// This function prints a single line of information that is
	// printed every so often (every UpdateScr) steps to inform the user
	// about the state of the simulation.
	void ArcusBase::infoLine() const
	{
		// prints a line of current state of system, monitors and energies

		// Information on the state of the protocol, i.e. step number, temperature,
		// etc...
		printf("%5d",Step );  

		// Information on the current data in the monitors
		mon.printCurData();

		// Information on the current energies
		ff->infoLine();
		
		// Newline (the previous calls don't print a newline character)
		printf("\n");
	}

	// this function should "match up" and label the columns produced by the
	// above function infoLine() and is usually called at the beginning of a run
	void ArcusBase::infoLineHeader() const
	{
		// prints the headers for the above function
		//printf("%5s", "NStep");

		// Headers of the monitors
		//mon.printHeader();

		// Headers of the forcefield components
		//ff->infoLineHeader();

		// New Line
		//printf("\n");
	}
}

