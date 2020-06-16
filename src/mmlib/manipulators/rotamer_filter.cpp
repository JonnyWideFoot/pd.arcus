#include "global.h"
#include <algorithm>
#include "rotamer_filter.h"

namespace Manipulator
{
	RotamerCumulativeProbabilityDensityFilter::RotamerCumulativeProbabilityDensityFilter( 
		const RotamerApplicatorBase& _app, 
		double _MinDensity )
		: RotamerFilterBase(_app),
		m_MinDensity( _MinDensity )
	{
		ASSERT(_MinDensity >= 0.0 && _MinDensity <= 1.0, ArgumentException, "Invalid fraction");
	}

	RotamerCumulativeProbabilityDensityFilter* RotamerCumulativeProbabilityDensityFilter::clone() const
	{
		return new RotamerCumulativeProbabilityDensityFilter(*this);
	}

	bool RotamerCumulativeProbabilityDensityFilter::passes( const RotamerLink& _ir, size_t _rotID )
	{
		double probability = _ir.probability( _rotID );
		if( probability <= 0.0 )
		{
			// Zero-probability or Undefined, so not active!
			return false;
		}

		// NOTE! This function does not and should not take notice of the active flag. 
		// We are just saying if this rotamer is or is not in the first X%
		size_t size = _ir.nRot();
		double probSum = 0.0;
		std::vector< std::pair< double, size_t > > probs;
		for( size_t i = 0; i < size; i++ )
		{
			double prob = _ir.probability( i );
			if( prob > 0.0 )
				probSum += prob; // -1.0 is the unassigned flag for a BBProb matrix			
			probs.push_back( std::pair< double, size_t >( prob, i ) );
		}
		std::sort( probs.rbegin(), probs.rend());
		double probCap = probSum * m_MinDensity;
		probSum = 0.0;
		for( size_t i = 0; i < size; i++ )
		{
			if( probs[i].second == _rotID )
				return true;
			probSum += probs[i].first;
			if( probSum > probCap )
				return false;
		}
		return false;
	}		


	RotamerDisulphideFilter::RotamerDisulphideFilter( const RotamerApplicatorBase& _app )
		: RotamerFilterBase(_app)
	{
	}

	RotamerDisulphideFilter* RotamerDisulphideFilter::clone() const
	{
		return new RotamerDisulphideFilter(*this);
	}

	bool RotamerDisulphideFilter::passes( const RotamerLink& _ir, size_t _rotID )
	{
		THROW( NotImplementedException, "RotamerDisulphideFilter not implemented");
	}


	RotamerCoordinationFilter::RotamerCoordinationFilter( const RotamerApplicatorBase& _app )
		: RotamerFilterBase(_app)
	{
	}

	RotamerCoordinationFilter* RotamerCoordinationFilter::clone() const
	{
		return new RotamerCoordinationFilter(*this);
	}

	bool RotamerCoordinationFilter::passes( const RotamerLink& _ir, size_t _rotID )
	{
		THROW( NotImplementedException, "RotamerCoordinationFilter not implemented");
	}
}

