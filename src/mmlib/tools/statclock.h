#ifndef __STATCLOCK_H
#define __STATCLOCK_H

#include <string>
#include <ctime>

//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API StatClock
{
	static const double TO_SECONDS;
	static const double TO_MILLI_SECONDS;

public:
	StatClock();
	StatClock( const std::string& name );
	void Reset();

	void Begin();
	void Stamp();
	void End();

	clock_t getAverage() const;
	clock_t getStdDev() const;
	clock_t getStdErr() const;
	clock_t getMax() const;
	clock_t getMin() const;

	bool HasData() const;
	std::string getName() const;

	void ReportSeconds( int _verbosity = 0 ) const;
	void ReportMilliSeconds( int _verbosity = 0 ) const;

	void setName( const std::string& _name ) { m_Name = _name; }

protected:
	// Internal Reporting function, called by public versions
	void Report( int _verbosity, double _TimeFactor, std::string _TimeName ) const;

	// Assertions
	void AssertBegun() const;
	void AssertNotBegun() const;
	void AssertHasData() const;

	// Member data
	std::string m_Name;
	bool m_Running;
	clock_t m_Last;
	std::vector<clock_t> m_Times;
};


//-------------------------------------------------
//
/// \brief Very silly 
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class PD_API TextProgressBar
{
public:
	TextProgressBar(int _Total, int _PrintWidth = 40);
	void Begin();
	void End();
	void next();
	void Reset(int _Total, int _PrintWidth);
	void Reset(int _Total);
	void Reset();

protected:
	void clearProgress();
	void DrawProgress();

	int m_Ticker;

	int m_PrintWidth;
	int m_Current;
	int m_Total;
};

#endif


