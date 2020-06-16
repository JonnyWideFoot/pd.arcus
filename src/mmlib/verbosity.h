#ifndef _VERBOSITY_H
#define _VERBOSITY_H

/// \brief An enumeration representing a PD-wide verbosity standard
/// \details Most objects have a member called 'Verbosity::Type OutputLevel;'
/// This defines the verbosity of that object. One then makes comparisons e.g.:
/// if( OutputLevel >= Verbosity::Normal ) Print("My wiiife!");
/// \author Jon Rea
class PD_API Verbosity
{
private:
	Verbosity();
public:
	enum Type
	{
		Silent = 0, ///< Just keep completely schtum!!
		Quiet = 2,  ///< Less output that normal. Just what I REALLY want to know!
		Normal = 4, ///< The kind of output you want to see every day
		Loud = 6,   ///< A bit more than normal. Mostly minor extra info.
		Scream = 8, ///< Debuggers paradise
		Eleven = 10 ///< Ever seen Spinal Tap? :-D God I'm whitty...
	};
};

#endif

