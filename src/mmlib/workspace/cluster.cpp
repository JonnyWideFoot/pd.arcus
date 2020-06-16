#include "global.h"
#include "workspace/cluster.h"


ClusterBase::ClusterBase():
	SnapShotLibrary()
{
}

ClusterBase::~ClusterBase()
{
}

size_t ClusterBase::clusterSize() const
{
	return cluster.size();
}

size_t ClusterBase::indexSize(const int indexNumber) const
{
	return cluster[indexNumber].index.size();
}

size_t ClusterBase::getCentre(const int indexNumber) const
{
	return cluster[indexNumber].centre;
}

size_t ClusterBase::getIndexMember(const int indexNumber, const int memberNumber) const
{
	return cluster[indexNumber].index[memberNumber];
}

void ClusterBase::reallocateCrmsArray()
{
	/// create a 2D array whose size in each dimension is initialised to the number of snapshots
	cRMSarray = TNT::Array2D <double> (dataSize(), dataSize());

	/// add the cRMS differences for each pair of SnapShots to the array
	for (int i = 0; i < dataSize(); i++)
	{
		for (int j = (i + 1); j < dataSize(); j++)
		{
			double cRMS = data[i].cRMSFrom(data[j]);
			cRMSarray[i][j] = cRMS;
			cRMSarray[j][i] = cRMS;
		}
		cRMSarray[i][i] = 0;
	}
}



ClusterSimple::ClusterSimple():
	ClusterBase()
{
}

ClusterSimple::~ClusterSimple()
{
}

void ClusterSimple::doClustering()
{
	///\brief do this each time the doClustering function is called to make sure the
	/// cRMS array is up to date
	reallocateCrmsArray();

	///\brief resize the taken array to make sure it fits the number of Snapshots
	/// and initialise everything to false
	taken.clear(); // resize() alone does not set all to false!
	taken.resize(dataSize(),false);

	///\brief do this each time the doClustering function is called to make sure the
	/// cluster vector is clean and ready to go
	cluster.clear();

	for (int i = 0; i < dataSize(); i++)
	{
		if (taken[i]) continue;

		///\brief make a new instance of an IndexGroup to hold each cluster while it is created
		/// then transfer its contents to the cluster vector at the end of the loop
		IndexGroup newCluster;

		///\brief if the SnapShot i hasn't already been assigned to an IndexGroup
		/// then make it the representative SnapShot of a new one
		newCluster.centre = i;

		for (int j = 0; j < dataSize(); j++)
		{
			if (taken[j]) continue;

			if (cRMSarray[i][j] < clusterCutoff)
			{
				newCluster.index.push_back(j);
				taken[j] = true;
			}
		}

		cluster.push_back( newCluster );
	}
}

