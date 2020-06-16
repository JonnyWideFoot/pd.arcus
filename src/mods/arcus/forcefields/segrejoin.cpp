#include "global.h"

#include "segrejoin.h"

using namespace Maths;

namespace Physics
{
	PeptideGroupRejoinForce::PeptideGroupRejoinForce(SegmentDef& _segdef, Verbosity::Type _verbosity)
		: SegmentDefUser(_segdef), // default constructor
		ForcefieldBase(_segdef.getWorkSpace()),
		m_ForceField(_segdef.getWorkSpace(),_verbosity)
	{
		Initialise();
		name = "PeptideGroupRejoinForce";
		m_CARestraintsOn = false;
	}

	bool PeptideGroupRejoinForce::Initialise( double strength, bool useCARestraints )
	{
		m_CARestraintsOn = useCARestraints;
		m_RestraintStrength = strength; // * Physics::PhysicsConst::kcal2JDivNa;
		setup(); // do all the internal jiggery-pokery
		return true;
	}

	inline void ForceIt( const WorkSpace& wspace, Physics::FF_Custom& parent, double strength, double idealDist, double gamma, size_t i, size_t j )
	{
		Physics::CustomForce cforceA( i, j );
		Physics::CustomForce cforceB( i, j );
		cforceA.setToVShaped(idealDist, strength, gamma,  1.0); // left arm
		cforceB.setToVShaped(idealDist, strength, gamma, -1.0); // right arm
		parent.addForce(cforceA);
		parent.addForce(cforceB);
	}

	void PeptideGroupRejoinForce::calcForces()
	{ 
		if( Active ) 
			m_ForceField.calcForces(); 
	}

	void PeptideGroupRejoinForce::enableHarmonicCARestraints( bool enabled )
	{
		if( enabled == true )
		{
			if( m_CARestraintsOn )
			{
				updateHarmonicCARestraints();
			}
			else
			{
				m_CARestraintsOn = true;
				setup();
			}
		}
		else	
		{
			m_CARestraintsOn = false;
			setup();
		}
	}

	const double DEFAULT_CA_WELL_WIDTH = 0.5;
	void PeptideGroupRejoinForce::addHarmonicCARestraints()
	{
		const WorkSpace& wspace = WorkSpaceOperatorBase::getWSpace();
		const SegmentDef& segDef = getSegDef();
		for( size_t i =	segDef.getStartResIndex(); i <= segDef.getEndResIndex(); i++ )
		{
			int caIndex = wspace.res[i].iCA;
			ASSERT( caIndex != -1, CodeException, "CA index lookup failed");
			const dvector& initialCAPos = wspace.atom[caIndex].pos();
			Physics::CustomForce cforceA( caIndex, initialCAPos );
			Physics::CustomForce cforceB( caIndex, initialCAPos );
			cforceA.setToVShaped(0.0 - DEFAULT_CA_WELL_WIDTH, m_RestraintStrength, 0.2,  1.0); // left arm
			cforceB.setToVShaped(0.0 + DEFAULT_CA_WELL_WIDTH, m_RestraintStrength, 0.2, -1.0); // right arm
			m_ForceField.addForce(cforceA);
			m_ForceField.addForce(cforceB);
		}
	}

	void PeptideGroupRejoinForce::updateHarmonicCARestraints()
	{
		ASSERT( m_CARestraintsOn == true, CodeException, "CA-Restraints are not configured, update impossible, call addHarmonicCARestraints()!!");
		const WorkSpace& wspace = WorkSpaceOperatorBase::getWSpace();
		const SegmentDef& segDef = getSegDef();
		size_t segLength = segDef.getNRes();
		size_t nRest = m_ForceField.size();
		ASSERT( segLength*2 < nRest, CodeException, "Insufficient CustomForces!?!" );
		for( size_t i =	segDef.getStartResIndex(); i <= segDef.getEndResIndex(); i++ )
		{
			int caIndex = wspace.res[i].iCA;
			ASSERT( caIndex != -1, CodeException, "CA index lookup failed");
			const dvector& updatedCAPos = wspace.atom[caIndex].pos();

			size_t a = nRest - 2 - ((segDef.getEndResIndex() - i)*2);
			size_t b = nRest - 1 - ((segDef.getEndResIndex() - i)*2);

			Physics::CustomForce& cforceA = m_ForceField.getForce(a);
			cforceA.setAbsPosTo( updatedCAPos );
			Physics::CustomForce& cforceB = m_ForceField.getForce(b);
			cforceB.setAbsPosTo( updatedCAPos );

			cforceA.setToVShaped(0.0 - DEFAULT_CA_WELL_WIDTH, m_RestraintStrength, 0.2,  1.0); // left arm
			cforceB.setToVShaped(0.0 + DEFAULT_CA_WELL_WIDTH, m_RestraintStrength, 0.2, -1.0); // right arm
		}
	}

	void PeptideGroupRejoinForce::setup()
	{
		const double gamma = 0.01; // rate of curve rounding....

		const WorkSpace& wspace = ForcefieldBase::getWSpace();
		const SegmentDef& segDef = getSegDef();
		m_ForceField.clear();		

		switch( segDef.getBreakType() )
		{
			case SegBreakEnd:
			{
				int resIndex = getSegDef().getEndResIndex();

				Physics::CustomForce cforceLC(
					wspace.res[ resIndex ].iC,
					wspace.res[ resIndex+1 ].iN );
				Physics::CustomForce cforceRC(
					wspace.res[ resIndex ].iC,
					wspace.res[ resIndex+1 ].iN );

				cforceLC.setToVShaped(OmegaGroupFilter::peptideDist_CN, m_RestraintStrength, gamma,  1.0); // left arm
				cforceRC.setToVShaped(OmegaGroupFilter::peptideDist_CN, m_RestraintStrength, gamma, -1.0); // right arm

				m_ForceField.addForce(cforceLC);
				m_ForceField.addForce(cforceRC);

				break;
			}
			case SegBreakStart:
			{
				int resIndex = getSegDef().getStartResIndex();

				Physics::CustomForce cforceLC(
					wspace.res[ resIndex ].iN,
					wspace.res[ resIndex-1 ].iC );
				Physics::CustomForce cforceRC(
					wspace.res[ resIndex ].iN,
					wspace.res[ resIndex-1 ].iC );

				cforceLC.setToVShaped(OmegaGroupFilter::peptideDist_CN, m_RestraintStrength, gamma,  1.0); // left arm
				cforceRC.setToVShaped(OmegaGroupFilter::peptideDist_CN, m_RestraintStrength, gamma, -1.0); // right arm

				m_ForceField.addForce(cforceLC);
				m_ForceField.addForce(cforceRC);

				break;
			}
			case SegBreakCentre:
			{
				int resIndex = getSegDef().getCentreResIndex();
				int resIndexP1 = resIndex+1;
				ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CN,  gamma, wspace.res[ resIndex ].iC, wspace.res[ resIndexP1 ].iN );
				ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_ON,  gamma, wspace.res[ resIndex ].iO, wspace.res[ resIndexP1 ].iN );
				ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CaN, gamma, wspace.res[ resIndex ].iCA, wspace.res[ resIndexP1 ].iN );
				ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CCa, gamma, wspace.res[ resIndex ].iC, wspace.res[ resIndexP1 ].iCA );				
				ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CaCa, gamma, wspace.res[ resIndex ].iCA, wspace.res[ resIndexP1 ].iCA );			
				if( wspace.res[ resIndexP1 ].iH != -1 )
				{
					// These will be ignored for proline, and will be INVALID for cis-residues
					ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_OH, gamma, wspace.res[ resIndex ].iO, wspace.res[ resIndexP1 ].iH );			
					ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CaH, gamma, wspace.res[ resIndex ].iCA, wspace.res[ resIndexP1 ].iH );			
					ForceIt( wspace, m_ForceField, m_RestraintStrength, OmegaGroupFilter::peptideDist_CH, gamma, wspace.res[ resIndex ].iC, wspace.res[ resIndexP1 ].iH );			
				}
				break;
			}
		default:
			{
				THROW( CodeException, "Unknown SegBreakType Encountered");
			}
		}

		if( m_CARestraintsOn )
		{
			addHarmonicCARestraints();
		}

		m_ForceField.setup(); // perform internal setup of the forcefield instance
	}
}


