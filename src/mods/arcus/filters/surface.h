#ifndef __SURFACE_H
#define __SURFACE_H

#include "filters/basicfilters.h"
#include "pickers/pickbase.h"
#include "tools/cloneholder.h"
#include "maths/maths_vector.h"
#include "maths/maths_matrix3x3.h"
#include "workspace/workspace.fwd.h"

class SegmentDefBase;

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class SurfacePicker : public PickBase
{
public:
	SurfacePicker();
	virtual ~SurfacePicker(); ///< You should always define a virtual destructor if you have virtual functions
	virtual SurfacePicker* clone() const;

	/// The box can be oriented purely based upon the anchors, however better placement 
	/// can be obtained by utilising the centre of geometry rather than the anchor centre
	/// as the frame of reference for the box orientation vector. It is recomented that this
	/// value remains at its default 'true'.
	bool COGMode; 

	virtual bool matches( const Particle& particle ) const;

	void calcBoundingBox( const SegmentDefBase& _seg );
	void drawBounds( IDrawProvider& _DrawHere );

private:
	const SegmentDefBase* m_Seg;
	CloneHolder<PickBase> m_AtomTypePicker; ///< We will pick this type of atom within the surface box. Defaults to CA atoms.

	// Particle transformation into the box below
	Maths::matrix3x3 sysRotation1, sysRotation2; ///< Define two rotations to place the box into position
	Maths::dvector systemDispacement; ///< Displacement to put anchor vector 1A onto the origin	

	// Internal box bounds
	double xMin; ///< Minimum coordinate bonds of the x axis
	double xMax; ///< Maximum coordinate bonds of the x axis
	double yMin; ///< Minimum coordinate bonds of the y axis
	double yMax; ///< Maximum coordinate bonds of the y axis
	double zMin; ///< Minimum coordinate bonds of the z axis
	double zMax; ///< Maximum coordinate bonds of the z axis
};






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka & Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class SegSurfaceFilter : public FilterBase
{
public:
	SegSurfaceFilter();
	SegSurfaceFilter( const SegmentDefBase& _seg );
	virtual ~SegSurfaceFilter();
	virtual SegSurfaceFilter* clone() const { return new SegSurfaceFilter(*this); }

	void setFilterTarget( const SegmentDefBase& _seg );

	virtual bool passes();
	virtual std::string reason();

private:
	void initCore();

	virtual void setInternalName();
	const SegmentDefBase* m_Seg;

	ClashFilter m_Clash;
};

#endif

