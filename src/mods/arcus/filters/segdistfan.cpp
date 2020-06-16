#include "global.h"
#include <istream>
#include "tools/io.h"
#include "tools/stringtool.h"
#include "workspace/segdef.h"
#include "workspace/snapshot.h"
#include "segdistfan.h"

SegmentDistanceFilter::Reach::Reach() 
	: min(DBL_MAX), 
	max(DBL_MIN), 
	CAPointer(NULL),
	AncPointer(NULL)
{
}

SegmentDistanceFilter::SegmentDistanceFilter( const std::string& _FileName ) 
	: FilterBase(),
	m_Mol(NULL),
	m_SegDef(NULL),
	m_FileName(_FileName)
{
	if( !IO::fileExists( _FileName ) )
	{
		throw IOException("SegmentDistanceFilter Could not find the specified input file");
	}
}

void SegmentDistanceFilter::setInternalName()
{
	name = "SegmentDistanceFan";
}

void SegmentDistanceFilter::initialise( const SegmentDefBase& _seg, bool _reversed )
{
	// Member assignment
	m_Mol = &_seg.getWorkSpace();
	m_SegDef = &_seg;
	m_Reversed = _reversed;

	// Setup
	anchorProperties(); // get the anchor separation and delta values
	obtainData(); // extract only a relevent subset of the data from a large input file
}

void SegmentDistanceFilter::anchorProperties()
{
	size_t nRes = m_Mol->res.size();
	ASSERT( m_SegDef->getStartResIndex() != 0 && nRes > m_SegDef->getEndResIndex(), 
		ProcedureException, "SegmentDistanceFilter: Cannot be applied to terminal segments" );
	
	int a1A = m_Mol->res[ m_SegDef->getStartResIndex()-1 ].iCA;
	int a1B = m_Mol->res[ m_SegDef->getStartResIndex()-1 ].iC;
	int a2A = m_Mol->res[ m_SegDef->getEndResIndex()+1 ].iCA;
	int a2B = m_Mol->res[ m_SegDef->getEndResIndex()+1 ].iN;
	if (a1A == -1 || a1B == -1 || a2A == -1 || a2B == -1) 
		THROW( ProcedureException, "SegmentDistanceFilter: Cannot find desired anchor atom." );

	SnapShotAtom* atom = m_Mol->cur.atom;

	m_AncSep = atom[a1A].p.dist(atom[a2A].p);
	double anchorSep2 = atom[a1B].p.dist(atom[a2B].p);
    m_AncDelta = m_AncSep - anchorSep2;

	m_OurSegLength = m_SegDef->getNRes();
	m_HalfOurSegLength = (m_OurSegLength + 1) / 2;

	if( m_Reversed )
	{
		// C-terminal
		size_t cnt = m_OurSegLength - m_HalfOurSegLength;
		m_Reach.clear();
		m_Reach.resize( cnt, Reach() );
		for( size_t i = 0; i < cnt; i++ )
		{
			// IMPORTANT!!!!!
			// NOTE, we are adding the FINAL residue of the loop into the STARTING m_Reach reach slot!!
			size_t at = m_SegDef->getEndResIndex() - i;
			size_t iCA = m_Mol->res[ at ].iCA;
			m_Reach[i].CAPointer = &atom[iCA].p;
			m_Reach[i].AncPointer = &atom[a1B].p;
		}
	}
	else
	{
		// N-terminal
		m_Reach.clear();
		m_Reach.resize( m_HalfOurSegLength, Reach() );
		for( size_t i = 0; i < m_HalfOurSegLength; i++ )
		{
			size_t at = m_SegDef->getStartResIndex() + i;
			size_t iCA = m_Mol->res[ at ].iCA;
			ASSERT( iCA != -1, CodeException, "CA lookup failure!" );
			m_Reach[i].CAPointer = &atom[iCA].p;
			m_Reach[i].AncPointer = &atom[a2B].p;
		}
	}
}

void SegmentDistanceFilter::obtainData()
{
	std::ifstream file(m_FileName.c_str(), std::ifstream::in);
	if( !file.is_open() ) 
		throw IOException( "SegmentDistFilter: Data file could not be opened: '" + m_FileName + "'!" );

	int restartCount = 0;
	const double AnchorSepToleranceDefault = 2.0; // Angstroms
	const double AnchorDeltaToleranceDefault = 0.2; // Angstroms
	const double AnchorSepToleranceDefaultGain = 0.2; // Angstroms
	const double AnchorDeltaToleranceDefaultGain = 0.02; // Angstroms
	const int MaxRestart = 15; // 2.0 + (0.2*MaxRestart) defines max sep...

	double AnchorSepTolerance = AnchorSepToleranceDefault;
	double AnchorDeltaTolerance = AnchorDeltaToleranceDefault;

RESTART_CUTOFF:
	if( restartCount++ >= MaxRestart )
	{
		THROW(ProcedureException,"SegmentDistanceFilter: Cutoff increase maxed out! Not enough data!");
	}
	int usingDataCount = 0;

	std::string line; // Temporary line container
	while(std::getline(file,line))
	{
		if( line.size() == 0 || line[0] == '#' )
		{
			continue; // dead line
		}

		std::vector<std::string> elements = chopstr( line, "," );
		size_t nElem = elements.size();

		int segLength = 0;
		double anchorSep = -1.0;
		double sepDelta = -1.0;

		ASSERT( nElem >= 4, ParseException, "SegmentDistFilter: Invalid line in input file");
		ASSERT( 0 == str2int(elements[0],segLength), ParseException, "SegmentDistFilter: Line is invalid, first element is not an integer");

		if( segLength != m_OurSegLength )
		{
			continue;
		}

		ASSERT( nElem == 3 + segLength, ParseException, "SegmentDistFilter: Invalid line in input file");
		ASSERT( 0 == str2double(elements[1],anchorSep), ParseException, "SegmentDistFilter: Line is invalid, anchorSep is not a float");
		ASSERT( 0 == str2double(elements[elements.size()-1],sepDelta), ParseException, "SegmentDistFilter: Line is invalid, sepDelta is not a float");

		double minSep = m_AncSep - AnchorSepTolerance;
		double maxSep = m_AncSep + AnchorSepTolerance;

		double minDelta = m_AncDelta - AnchorDeltaTolerance;
		double maxDelta = m_AncDelta + AnchorDeltaTolerance;

		if( minSep > anchorSep || anchorSep > maxSep )
		{
			continue;
		}

		if( minDelta > sepDelta || sepDelta > maxDelta )
		{
			continue;
		}

		// Monitor how much data we have
		usingDataCount++;

		// We have found a PDB derived loop data set, 
		// which matches this loops length, anchor separation, and anchor orientation!
		// Use this data below ...
		size_t to = (size_t) segLength;
		for( size_t i = 0; i < to; i++ )
		{
			double dist;
			ASSERT( 0 == str2double(elements[2+i],dist), ParseException, "SegmentDistFilter: Line is invalid, element is not a float");

			size_t at = i;
			if( i >= m_HalfOurSegLength )
			{
				at = segLength - i - 1;
			}

			if( at != m_Reach.size() ) // happens when we have an odd number of residues in the loop...
			{
				Reach& r = m_Reach[at];
				r.max = std::max( dist, r.max );
				r.min = std::min( dist, r.min );
			}
		}
	}

	const int REQUIRED_DATA = 15;
	if( usingDataCount < REQUIRED_DATA )
	{
		Printf("SegmentDistanceFilter: Cutoffs are too tight! Adapting...\n");

		// back to the start of the file!!		
		file.clear(); // first clear the eof and fail bits
		file.seekg(0, std::ios::beg); // then seek...

		// Increase cutoffs
		AnchorSepTolerance += AnchorSepToleranceDefaultGain;
		AnchorDeltaTolerance += AnchorDeltaToleranceDefaultGain;

		goto RESTART_CUTOFF;
	}
	Printf("SegmentDistanceFilter: Initialised OK!\n  AnchorSepTollerance: %6.3lf (%6.3lf)\n  AnchorSepDeltaTollerance: %6.3lf (%6.3lf)\n")
		(AnchorSepTolerance)(AnchorSepToleranceDefault)(AnchorDeltaTolerance)(AnchorDeltaToleranceDefault);

	double scaleDataBy = 0.05; // 5%
	for( size_t i = 0; i < m_Reach.size(); i++ )
	{
		ASSERT( m_Reach[i].max != DBL_MIN && m_Reach[i].min != DBL_MAX,
			CodeException, "If we have data, which should have already been asserted, this should not happen!");
		ASSERT( m_Reach[i].max >= 0.0 && m_Reach[i].min >= 0.0,
			CodeException, "Ewww! Anomaly!!");
		// Store as SQUARED values!!
		m_Reach[i].max += 0.01; // 0.01 for rounding errors in the self comparison...
		m_Reach[i].min -= 0.01; // 0.01 for rounding errors in the self comparison...
		m_Reach[i].max += m_Reach[i].min * scaleDataBy; // Allow for experimental wobble
		m_Reach[i].min -= m_Reach[i].min * scaleDataBy;
		m_Reach[i].max *= m_Reach[i].max; 
		m_Reach[i].min *= m_Reach[i].min;
	}

	file.close();
}

void SegmentDistanceFilter::draw( IDrawProvider& _DrawMe )
{
	ASSERT( m_Mol != NULL, UninitialisedException, "SegmentDistanceFilter used prior to init!");
	DrawingVector* dv;
	for( size_t i = 0; i < m_Reach.size(); i++ )
	{
		if( NULL != (dv = _DrawMe.request()) )
		{
			double sqrDist = m_Reach[i].sqrdist();
			if( sqrDist < m_Reach[i].min || sqrDist > m_Reach[i].max )
			{
				dv->colourCode = Colour::Red;
			}
			else
			{
				dv->colourCode = Colour::Green;
			}
			dv->v1.setTo( *m_Reach[i].CAPointer );
			dv->v2.setTo( *m_Reach[i].AncPointer );
		}
		else
		{
			return;
		}
	}
}

bool SegmentDistanceFilter::passes()
{
	ASSERT( m_Mol != NULL, UninitialisedException, "SegmentDistanceFilter used prior to init!");
	for( size_t i = 0; i < m_Reach.size(); i++ )
	{
		double sqrDist = m_Reach[i].sqrdist();

		//double dist = m_Reach[i].dist();
		//double bobMin = sqrt(m_Reach[i].min);
		//double bobMax = sqrt(m_Reach[i].max);

		if( sqrDist < m_Reach[i].min || sqrDist > m_Reach[i].max )
		{
			return false;
		}
	}
	return true;
}

std::string SegmentDistanceFilter::reason()
{
	ASSERT( m_Mol != NULL, UninitialisedException, "SegmentDistanceFilter used prior to init!");
	for( size_t i = 0; i < m_Reach.size(); i++ )
	{
		double sqrDist = m_Reach[i].sqrdist();
		if( sqrDist < m_Reach[i].min || sqrDist > m_Reach[i].max )
		{
			size_t res = m_Reversed ? m_SegDef->getEndResIndex() - i : m_SegDef->getStartResIndex() + i;
			StringBuilder sb;
			sb.setFormat("Residue %d CA-atom is out of its bounds. Untrue: %6.3f <= %6.3f <= %6.3f")
				(res)(sqrt(m_Reach[i].min))(sqrt(sqrDist))(sqrt(m_Reach[i].max));
			return sb.toString();
		}
	}
	return "No over-reach, all fine.";
}
