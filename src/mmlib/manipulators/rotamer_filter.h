#ifndef __ROTAMER_FILTER
#define __ROTAMER_FILTER

#include "rotamer_applicatorbase.h" // base class

namespace Manipulator
{
	/// \brief A filter that ensures that only the most likely states are kept within the library.
	/// \details SCWRL 3.0 (original paper, page 2005, column 1, step 2, states that "rotamers are read in from
	/// highest to lowest probability until the cumulative density reaches at least 90%." I take that to mean that
	/// the <10% of low probability rotamers are removed. This class acomplishes that in any given library by
	/// filtering these improbable states.
	class PD_API RotamerCumulativeProbabilityDensityFilter : public RotamerFilterBase
	{
	public:
		RotamerCumulativeProbabilityDensityFilter( const RotamerApplicatorBase& _app, double _MinDensity );
		virtual ~RotamerCumulativeProbabilityDensityFilter(){}
		virtual RotamerCumulativeProbabilityDensityFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
		double getMinDensity() const { return m_MinDensity; }
	private:
		double m_MinDensity;
	};

	/// Only allow CYS rotamers which are compatible with forming a disulphide bond
	class PD_API RotamerDisulphideFilter : public RotamerFilterBase
	{
	public:
		RotamerDisulphideFilter( const RotamerApplicatorBase& _app );
		virtual ~RotamerDisulphideFilter(){}
		virtual RotamerDisulphideFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
	};

	/// Only allow rotamers that are compatible with the coordination of a given location in space.
	class PD_API RotamerCoordinationFilter : public RotamerFilterBase
	{
	public:
		RotamerCoordinationFilter( const RotamerApplicatorBase& _app );
		virtual ~RotamerCoordinationFilter(){}
		virtual RotamerCoordinationFilter* clone() const;
		virtual bool passes( const RotamerLink& _ir, size_t _rotID );
	};
}

#endif

