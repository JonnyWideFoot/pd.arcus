#ifndef __CLUSTER_H
#define __CLUSTER_H

#include <vector>

#include "maths/tntjama/tnt_array2d.h"
#include "maths/tntjama/tnt_array2d_utils.h"

#include "workspace/snapshot.h"

struct IndexGroup
{
	std::vector <size_t> index;
	size_t centre;
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
class ClusterBase: public SnapShotLibrary
{
public:

	/// constructor
	ClusterBase();

	/// destructor
	virtual ~ClusterBase();

	/// pure virtual function that houses the clustering algorithm
	virtual void doClustering() = 0;

	/// get the size of the std::vector cluster, a protected member variable
	size_t clusterSize() const;

	/// get the size of the std::vector index found at cluster[indexNumber]
	size_t indexSize(const int indexNumber) const;

	///\brief get the number of the representative SnapShot (centre) of the std::vector index
	/// found at cluster[indexNumber
	size_t getCentre(const int indexNumber) const;

	/// get the number of the SnapShot found at cluster[indexNumber].index[memberNumber]
	size_t getIndexMember(const int indexNumber, const int memberNumber) const;

	/// allows you to set the cRMS cutoff for clustering
	double clusterCutoff;

protected:
	std::vector <IndexGroup> cluster;

	/// 2D Array containing cRMS values for each pair of SnapShots
	TNT::Array2D <double> cRMSarray;

	///\brief function to initialise a 2D array to the number of SnapShots and
	/// add to the array cRMS differences between each pair of SnapShots
	virtual void reallocateCrmsArray();

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
class ClusterSimple: public ClusterBase
{
public:
	/// constructor
	ClusterSimple();

	/// destructor
	virtual ~ClusterSimple();

	/// NB. not pure virtual
	virtual void doClustering();

protected:

	///\brief an array of bools allows the clustering function to keep track of
	/// which SnapShots have already been assigned to a cluster
	std::vector <bool> taken;
};

#endif
