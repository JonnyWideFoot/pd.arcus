#ifndef __SEG_DIST_FAN_H
#define __SEG_DIST_FAN_H

#include <string>
#include "filters/filterbase.h"
#include "workspace/workspace.fwd.h"

// Forward declarations
class SegmentDefBase;

//-------------------------------------------------
/// \brief  BRIEF DESCRIPTION
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// \author Jon Rea 
/// \todo STATE OF DEVELOPMENT
/// \bug BUGS?
class SegmentDistanceFilter : public FilterBase
{
public:
	SegmentDistanceFilter( const std::string& _FileName );
	virtual ~SegmentDistanceFilter(){}
	virtual SegmentDistanceFilter* clone() const { return new SegmentDistanceFilter(*this); }

	void initialise( const SegmentDefBase& _seg, bool _reversed );

	virtual bool passes();
	virtual std::string reason();

	void draw( IDrawProvider& _DrawMe );

protected:
	virtual void setInternalName();
	const SegmentDefBase* m_SegDef;
	WorkSpace* m_Mol;
	bool m_Reversed;

private:
	void anchorProperties();
	void obtainData();
	std::string m_FileName;
	double m_AncSep;
	double m_AncDelta;
	size_t m_OurSegLength;
	size_t m_HalfOurSegLength;

	struct Reach
	{
		Reach();
		double min;
		double max;
		const Maths::dvector* CAPointer;
		const Maths::dvector* AncPointer;
		inline double sqrdist() const { return CAPointer->sqrdist( *AncPointer ); }
		inline double dist() const { return CAPointer->dist( *AncPointer ); }
	};
	std::vector< Reach > m_Reach;
};

#endif

