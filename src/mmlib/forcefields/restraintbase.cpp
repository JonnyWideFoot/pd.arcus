#include "global.h"

#include "workspace/workspace.h"
#include "workspace/neighbourlist.h"
#include "pickers/basicpickers.h"

#include "forcefields/forcefield.h"
#include "forcefields/restraintbase.h"

// namespace includes
using namespace Maths;

namespace Physics
{

	RestraintForcefieldBase::RestraintForcefieldBase(
		WorkSpace &newwspace, 
		double _k
	):
		ForcefieldBase( newwspace )
	{
		k = _k;
		m_RestraintStruc_set = false;        // restraint structre is not set by default
	}


	void RestraintForcefieldBase::setRestraintStructure(const SnapShot &psp)
	{
		m_RestraintStruc     = psp;
		m_RestraintStruc_set = true; // flag that a psp structure has been deposited
		needsetup            = true; // flag that we will need to be re-setup
	}

	void RestraintForcefieldBase::setup()
	{		
		SnapShot savesystem;

		if(m_RestraintStruc_set)
		{
			// preserve current workspace configuration
			savesystem = getWSpace().save();
		}

		// if the user has set a particular structure to restrain, load that
		if(m_RestraintStruc_set)
		{
			getWSpace().load_forced(m_RestraintStruc);
		}

		// remember current system setup
		saveCurrentAtomPositions();

		// if the user had set a particular structure to restrain to
		// the restore the original system position and remove the flag that
		// the user has set something
		if(m_RestraintStruc_set)
		{
			getWSpace().load(savesystem);
			m_RestraintStruc_set = false;
		}

	}

	void RestraintForcefieldBase::saveCurrentAtomPositions()
	{
		// save the atom positions.
		printf("Setting Restraint position for: %s \n",name.c_str());
		getWSpace().outtra.append();

	}










	AtomRestraintForcefieldBase::AtomRestraintForcefieldBase(
		WorkSpace &newwspace, 
		double _k
	):
		RestraintForcefieldBase( newwspace, _k )
	{
		m_Picker = PickAllParticles();   // restraint applies to all atoms by default
	}

	void AtomRestraintForcefieldBase::setSelection(const PickBase &_Picker)
	{ 
		// set the new user's custom picker of choice
		m_Picker = _Picker; // copy over a clone of the supplied picker
		needsetup = true;
	}


	void AtomRestraintForcefieldBase::setup(){
		// just use base class's etup() function
		RestraintForcefieldBase::setup();
	}

	void AtomRestraintForcefieldBase::saveCurrentAtomPositions()
	{
		RestraintForcefieldBase::saveCurrentAtomPositions();

		// Create the selection using the current Picker
		m_Selection.setPicking(getWSpace(),m_Picker.data());
		m_Selection.store(); // ensure we store the current positions!
	}


	void  AtomRestraintForcefieldBase::detail() 
	{
		setup();
		info();
		// print which atoms will be positionally restrained
		printf(" Details of atoms restrained: \n\n" ); 
		getWSpace().detail( m_Selection.getPicker() );
		printf("\n\n" ); 
	}



} // namespace Physics


