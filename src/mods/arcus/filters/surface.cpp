#include "global.h"
#include "workspace/workspace.h"
#include "workspace/segdef.h"
#include "surface.h"

using namespace Maths; // Shortens some of the commands below

SurfacePicker::SurfacePicker() : m_Seg( NULL ), COGMode(true)
{
	m_AtomTypePicker = PickCAAtoms();
}

SurfacePicker::~SurfacePicker()
{
}

SurfacePicker* SurfacePicker::clone() const
{
	return new SurfacePicker(*this);
}

void SurfacePicker::calcBoundingBox( const SegmentDefBase& _seg )
{
	// Check resnl.
	const int checkDecimalPlaces = 8;

	// Proxies
	m_Seg = & _seg;
	WorkSpace& wspace = _seg.getWorkSpace();
	int startResIndex = _seg.getStartResIndex();
	int segLength = _seg.getNRes();
	int lastLoopResIndexPlus1 = startResIndex + segLength;

	ASSERT( startResIndex != 0 && lastLoopResIndexPlus1 < (int)wspace.nResidues(), ArgumentException, "A surface box cannot be calculated for a terminal segment at the moment!");

	// -------------------------------------------------------------------
	// PHASE 1: Get the atom positions and the translation onto the origin
	// We will use the position the midpoint between the two anchor points
	// as the new origin: "vAncMA". The translation will be held in the
	// local variable "systemDispacement".
	// -------------------------------------------------------------------

	// get the anchor atom vectors
	// these are defined as the N and C of the peptide groups joing the loop to the protein
	dvector vAnc1A, vAnc1B, vAnc2A, vAnc2B;
	vAnc1A.setTo( wspace.atomxyz( wspace.res[ startResIndex - 1 ].iC ));
	vAnc1B.setTo( wspace.atomxyz( wspace.res[ startResIndex ].iN ));
	vAnc2A.setTo( wspace.atomxyz( wspace.res[ lastLoopResIndexPlus1 ].iN ));
	vAnc2B.setTo( wspace.atomxyz( wspace.res[ lastLoopResIndexPlus1 - 1 ].iC ));

	// calculate the centre dvector between the two anchors
	dvector vFrom1AToMidPoint, vAncMB;
	vFrom1AToMidPoint.setTo(vAnc2A);
	vFrom1AToMidPoint.sub(vAnc1A);
	vFrom1AToMidPoint.div(2.0);

	// calculate the displacement to position the entire system on the new origin, this is centred on vAncMA
	systemDispacement.setTo(vAnc1A);
	systemDispacement.add(vFrom1AToMidPoint);

	vAnc1A.sub(systemDispacement);
	vAnc1B.sub(systemDispacement);
	vAnc2A.sub(systemDispacement);
	vAnc2B.sub(systemDispacement);
	// vAncMA.sub(systemDispacement); vAncMA is now == 0,0,0

	// Calculate the orientation vector vAncMB which points away from the box's cartesian-system origin
	if( COGMode )
	{
		// Centre of geometry box orientation mode
		dvector cog = wspace.getCentreOfGeometry();
		vAncMB.setTo(systemDispacement);	
		vAncMB.sub(cog);
		vAncMB.unify();
	}
	else
	{
		// Anchor box orientation mode
		vAncMB.setTo(vAnc1B);
		vAncMB.add(vFrom1AToMidPoint);
		vAncMB.sub(systemDispacement);
	}

	// -------------------------------------------------------------------
	// PHASE 2: Calculate the first rotation to place the dvector pointing
	// from atom A of anchor 1 to atom A of anchor 2 along the Y-axis
	// -------------------------------------------------------------------

	dvector vUnitX( 1.0, 0.0, 0.0 ); // xAxis unit dvector
	dvector vUnitY( 0.0, 1.0, 0.0 ); // yAxis unit dvector

	dvector vACrossU; // calculate the unit dvector of the midpoint of the anchors
	vACrossU.setTo(vAnc2A);
	///vACrossU.sub(vAncMA); ignore ... effectively mathematically sub(0,0,0)
	vACrossU.unify();

	if( !SigFigEquality( vACrossU.sqrdist(vUnitY), 0.0, checkDecimalPlaces ) )
	{
		double calcAngle = acos( dotProduct(vACrossU, vUnitY) );
		dvector vNormal; // normal to the plane formed by the above two vectors
		vNormal.crossProduct(vACrossU, vUnitY);
		sysRotation1.setToAxisRot(vNormal, calcAngle);// get the roatation
	}
	else
	{
		sysRotation1.setToIdentity();
	}

	//vAnc1A.mulmat( &sysRotation1); // Required for DEBUG only
	//vAnc1B.mulmat( &sysRotation1); // Required for DEBUG only
	//vAnc2A.mulmat( &sysRotation1); // Required for DEBUG only
	//vAnc2B.mulmat( &sysRotation1); // Required for DEBUG only
	vAncMB.mulmat(sysRotation1);

	// -------------------------------------------------------------------
	// PHASE 3: Calculate the second rotation to align the vAncMA to vAncMB
	// dvector with the plane made my the Y-axis.
	// -------------------------------------------------------------------

	dvector vYAxizToMBU;
	vYAxizToMBU.x = vAncMB.x;
	vYAxizToMBU.y = 0;
	vYAxizToMBU.z = vAncMB.z;
	vYAxizToMBU.unify();

	// Prevent maths exceptions by testing the square distance between unit vecors
	double dTemp = vUnitX.sqrdist(vYAxizToMBU);
	if( !SigFigEquality( dTemp, 0.0, checkDecimalPlaces ) )
	{
		double dotP = dotProduct(vUnitX, vYAxizToMBU);
		double calcAngle = acos( dotP );
		sysRotation2.setToYrot(calcAngle);
		vYAxizToMBU.mulmat( sysRotation2 );

		// verify that we have rotated in the correct direction
		// if we have rotated around the Y-axis in the wrong direction,
		// then set the rotation matrix to rotate in the other direction		
		if( !SigFigEquality( vUnitX.sqrdist(vYAxizToMBU), 0.0, checkDecimalPlaces ) )
		{
			sysRotation2.setToYrot( -calcAngle );

			// reverse the previous rotation twice, to rotate the other way
			vYAxizToMBU.mulmat(sysRotation2);
			vYAxizToMBU.mulmat(sysRotation2);

			// CHECK!!
			if( !SigFigEquality( vUnitX.sqrdist(vYAxizToMBU), 0.0, checkDecimalPlaces ) )
			{
				THROW(CodeException,"Maths vector error in calcBoundingBox()");
			}
		}
	}
	else
	{
		// It already is aligned along the x axis - not very likely, but possible
		sysRotation2.setToIdentity();
	}

	// -------------------------------------------------------------------
	// PHASE 4: Calculate the rerquired dimensions of the box that we are
	// drawing around the anchor points. Any CA atoms within this box will
	// be deemed to be surface atoms and included in the list ...
	// -------------------------------------------------------------------

	// calculate a very simplified "maximum loop extension" box around the anchors

	double distA1A2 = vAnc1A.dist(vAnc2A);
	// 'maxExtension' longest possible distance in extended conformation ( obviously too long for the box )
	// 3.8 is the CA-CA distance
	double maxExtension = ((double)segLength) * 3.8;
	if( maxExtension < distA1A2 )
	{
		throw ProcedureException("Surface analysis failure; The loop is too short to span the gap between the anchors!\n");
	}
	maxExtension *= 0.80; // Empiricle shrink factor

	// calculated from solving a right angle triangle pointing away from the Y-axis
	double maxReach = 0.5 * (maxExtension - ((distA1A2*distA1A2) / maxExtension));
	xMax = zMax = maxReach;
	xMin = zMin = -maxReach;
	xMin *= 0.3; // empiracle shrink factor of xMin only -  shallower cut into the protein core

	// the A1A2 dvector lies along the Y-axis by definition above - gives us less of the core atoms
	// This does take 'distA1A2' into account, but the term mathematically cancels out.
	yMax = maxExtension / 2.0;
	yMin = -yMax;
}

void SurfacePicker::drawBounds( IDrawProvider& _DrawHere )
{
	ASSERT( m_Seg != NULL, CodeException, "calcBoundingBox needs be be called before this instance is initialised");

	DrawingVector *vect = NULL;
	matrix3x3 tempRot;

	dvector b1(xMin, yMin, zMin);
	dvector b2(xMin, yMax, zMin);
	dvector b3(xMin, yMin, zMax);
	dvector b4(xMin, yMax, zMax);
	dvector t1(xMax, yMin, zMin);
	dvector t2(xMax, yMax, zMin);
	dvector t3(xMax, yMin, zMax);
	dvector t4(xMax, yMax, zMax);

	tempRot.setTo( sysRotation2 );
	tempRot.transpose();
	b1.mulmat(tempRot);
	b2.mulmat(tempRot);
	b3.mulmat(tempRot);
	b4.mulmat(tempRot);
	t1.mulmat(tempRot);
	t2.mulmat(tempRot);
	t3.mulmat(tempRot);
	t4.mulmat(tempRot);

	tempRot.setTo( sysRotation1 );
	tempRot.transpose();
	b1.mulmat(tempRot);
	b2.mulmat(tempRot);
	b3.mulmat(tempRot);
	b4.mulmat(tempRot);
	t1.mulmat(tempRot);
	t2.mulmat(tempRot);
	t3.mulmat(tempRot);
	t4.mulmat(tempRot);

	b1.add( systemDispacement );
	b2.add( systemDispacement );
	b3.add( systemDispacement );
	b4.add( systemDispacement );
	t1.add( systemDispacement );
	t2.add( systemDispacement );
	t3.add( systemDispacement );
	t4.add( systemDispacement );

	if( vect = _DrawHere.request() ) 
	{ vect->colourCode = Colour::Blue; vect->v1.setTo( fvector(b1) ); vect->v2.setTo( fvector(b2) ); }
	if( vect = _DrawHere.request() ) 
	{ vect->colourCode = Colour::Blue; vect->v1.setTo( b1 ); vect->v2.setTo( b3 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Blue; vect->v1.setTo( b4 ); vect->v2.setTo( b2 ); }				
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Blue; vect->v1.setTo( b4 ); vect->v2.setTo( b3 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Yellow; vect->v1.setTo( t1 ); vect->v2.setTo( t2 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Yellow; vect->v1.setTo( t1 ); vect->v2.setTo( t3 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Yellow; vect->v1.setTo( t4 ); vect->v2.setTo( t2 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Yellow; vect->v1.setTo( t4 ); vect->v2.setTo( t3 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Orange; vect->v1.setTo( b1 ); vect->v2.setTo( t1 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Orange; vect->v1.setTo( b2 ); vect->v2.setTo( t2 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Orange; vect->v1.setTo( b3 ); vect->v2.setTo( t3 ); }
	if( vect = _DrawHere.request() )
	{ vect->colourCode = Colour::Orange; vect->v1.setTo( b4 ); vect->v2.setTo( t4 ); }
}

bool SurfacePicker::matches( const Particle& particle ) const
{
	ASSERT( m_Seg != NULL, CodeException, "calcBoundingBox() needs be be called before this instance is initialised");

	// ---------------------------------------------------------------------------
	// PHASE 5: Check through all the m_AtomTypePicker atoms that are not involved
	// in the loop itself to see if they should be deemed "surface atoms"
	// ---------------------------------------------------------------------------

	// Proxies
	WorkSpace& wspace = m_Seg->getWorkSpace();
	const PickBase& picker = m_AtomTypePicker.data();

	if( particle.i >= (int)m_Seg->getStartAtomIndex() && particle.i <= (int)m_Seg->getEndAtomIndex() )
	{
		// It is within the segment and therefore cant be part of the surface
		return false;
	}
	
	if( !picker.matches( wspace.atom[particle.i] ) ) 
	{
		// only include atoms also included by the type picker
		return false; 
	}

	dvector checkVector( wspace.atomxyz(particle.i) );
	checkVector.sub( systemDispacement );
	checkVector.mulmat( sysRotation1 );
	checkVector.mulmat( sysRotation2 );	

	// Return if that the atom is inside the surface box
	return ( 
		checkVector.x > xMin && checkVector.x < xMax &&
		checkVector.y > yMin && checkVector.y < yMax &&
		checkVector.z > zMin && checkVector.z < zMax 
		);
}

SegSurfaceFilter::SegSurfaceFilter()
	: FilterBase(), 
	m_Seg(NULL)
{
}

SegSurfaceFilter::SegSurfaceFilter( const SegmentDefBase& _seg ) 
	: FilterBase(), 
	m_Seg(NULL)
{
	setFilterTarget( _seg);
}

SegSurfaceFilter::~SegSurfaceFilter()
{
}

void SegSurfaceFilter::setFilterTarget( const SegmentDefBase& _seg )
{
	m_Seg = &_seg;
	ASSERT( m_Seg != NULL, CodeException, "Internal NULL reference found in SegSurfaceFilter during setFilterTarget()");
	initCore();
}

void SegSurfaceFilter::initCore()
{
	// We are going to setup m_Clash to do the steric filtering
	ASSERT( m_Seg != NULL, CodeException, "Internal NULL reference found in SegSurfaceFilter. Has setFilterTarget() been called?");

	// Pick for the CA atoms in the loop
	m_Clash.setForPicker( Pick_AND( PickResidueRange(*m_Seg), PickCAAtoms() ) );

	// Pick against the CA atoms on the protein surface
	SurfacePicker picker;
	picker.calcBoundingBox( *m_Seg );
	m_Clash.setAgainstPicker( picker );

	// Finally set the molecule that the surface filter points to
	m_Clash.setMolecule( m_Seg->getWorkSpace() );
}

bool SegSurfaceFilter::passes()
{
	ASSERT( m_Seg != NULL, CodeException, "Internal NULL reference found in SegSurfaceFilter. Has setFilterTarget() been called?");
	return m_Clash.passes();
}

std::string SegSurfaceFilter::reason()
{
	ASSERT( m_Seg != NULL, CodeException, "Internal NULL reference found in SegSurfaceFilter. Has setFilterTarget() been called?");
	return m_Clash.reason();
}

void SegSurfaceFilter::setInternalName()
{
	name = "SegSurfaceFilter";
}

