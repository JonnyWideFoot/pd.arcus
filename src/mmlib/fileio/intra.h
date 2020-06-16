#ifndef __INPUTTRA_H
#define __INPUTTRA_H

#include <string>
#include "object.h"
#include "workspace/snapshot.fwd.h" 

/// \brief Basic base clas for trajectories - this is geared around reading in
/// typical MM trajectories which contain only coordinate data, not topology data (unlike the UOB tra format)
///
/// \details
/// It interfaces with SnapShot since that is the internal representation of such a data set
/// The basic version is geared only towards sequential reading, not random access. 
/// A derived base class also allows random access reading.
///
/// \author Mike Tyka
///

class PD_API InputTrajectory : public Object
{
public:
	InputTrajectory(const std::string& _filename) { filename = _filename; }
	virtual ~InputTrajectory () {}
	virtual InputTrajectory* clone() const = 0;

	/// \brief Reads a structure and returns a SnapShot with the coordinates (and the box geometry)
	/// \return Returns true if end of file was reached during read, false otherwise.
	virtual bool readNext( SnapShot &ss ) = 0;

	/// \brief Skips an entry in the file 
	/// \return Returns true if end of file was reached during read, false otherwise.
	virtual bool skip() = 0;

	/// \brief Returns true if file handle has reached end of file. 
	/// \return Returns true if end of file was reached during read, false otherwise.
	virtual bool isEndOfFile() const = 0;

	/// \brief Go back to the start of the tra
	virtual void reset() = 0;

protected:

	std::string filename;
};


class PD_API InputTrajectory_RandomAccess: public InputTrajectory
{
public:
	InputTrajectory_RandomAccess(const std::string& _filename) : InputTrajectory(_filename) {}
	virtual ~InputTrajectory_RandomAccess () {}
	virtual InputTrajectory_RandomAccess* clone() const = 0;

	/// \brief Read an arbitrary entry fromt he trajecotry 
	/// \return Returns true if end of file was reached during read, false if read was successful 
	virtual void readRandomAccess( SnapShot &ss, size_t entry ) = 0;

	/// \brief works out how many entries are in the trajectory
	virtual size_t nEntries() const = 0;

protected:
};

#endif


