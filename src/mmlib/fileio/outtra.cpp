#include "global.h"
#include "outtra.h"

#include "system/system.h"
#include "system/molecule.h"
#include "workspace/workspace.h"
#include "workspace/componentbase.h"

namespace IO
{

	// ------------------
	//  OutputTrajectory
	// ------------------

	OutputTrajectory::OutputTrajectory(  
		WorkSpace &_wspace 
	) :
		WorkSpaceOperatorBase( _wspace ),
		created(false)
	{

	}

	// ----------------------
	//  OutputTrajectoryFile
	// ----------------------

	OutputTrajectoryFile::OutputTrajectoryFile(
		const std::string& _filestem, 
		WorkSpace & _wspace
	) :
		filestem(_filestem), 
		OutputTrajectory( _wspace )
	{

	}

	// ---------------------------
	//  OutputTrajectoryContainer
	// ---------------------------


	OutputTrajectoryContainer::OutputTrajectoryContainer( 
		WorkSpace & _wspace 
	) :
		OutputTrajectory( _wspace )
	{

	}


	void OutputTrajectoryContainer::prepare()
	{
		if( size() > 0 ) 
		{
			getWSpace().calcCRMS_HeavyAtom();
		}
	}


	int OutputTrajectoryContainer::create()
	{
		prepare();
		for(unsigned i=0;i<size();i++) 
			element(i).create();
		return 0;
	}


	int OutputTrajectoryContainer::append()
	{
		prepare();
		for(unsigned i=0;i<size();i++) 
		{
			element(i).append();
		}
		return 0;
	}







	OutputFile::OutputFile(const std::string& _filestem)
	{
		filestem = _filestem;
	}

	void OutputFile::save( System &_system){
		WorkSpace wspace( _system );
		save( wspace );
	}

	void OutputFile::save( Molecule &_molecule){
		System mysystem( _molecule.ffps() );
		mysystem.add( _molecule );
		save( mysystem );
	}


} // namespace IO


