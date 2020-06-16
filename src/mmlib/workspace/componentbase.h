#ifndef __COMPONENT_BASE_H
#define __COMPONENT_BASE_H

#include "workspace/workspace.fwd.h"

//-------------------------------------------------
/// \brief  Base class for everything that needs to know and is permanently connected to WorkSpace
/// \details Operators are "things that work on workspace and live outside it"
/// \author Mike Tyka
class PD_API WorkSpaceOperatorBase
{
	friend class WorkSpace;
public:		
	WorkSpaceOperatorBase(WorkSpace &_wspace): m_WSpace( &_wspace){}

	WorkSpace       &getWSpace() const { return *m_WSpace; }

private:
	/// internal pointer to keep the workspace pointer. Once created this cannot be changed.
	WorkSpace *m_WSpace;
};



//-------------------------------------------------
/// \brief  Base class for everything that has a dynamic pointer to a Workspace and can be reinitialised.
/// \details Components are "things that work on workspace but live inside it (like a plug in)"
/// \author Mike Tyka
class PD_API WorkSpaceComponentBase
{
	friend class WorkSpace;
public:
	/// initialise into unconnected state
	WorkSpaceComponentBase(): wspace(NULL){}

	/// initialise into connected state
	WorkSpaceComponentBase(WorkSpace *_wspace): wspace(_wspace){}

	const WorkSpace *getWSpace() const { return wspace; }

protected:
	/// This function should not be overloaded by deriving classes
	/// it allows WorkSpace to call reinit(...) in derived classes
	/// without having to define a "friend class WorkSpace" in each daughter class
	void reinit_base( WorkSpace* _wspace )
	{
		reinit(_wspace);
	}

	/// Reassign to a new workspace - Derived classes
	/// can (and often must) implement more specific versions of this.
	virtual void reinit( WorkSpace* _wspace )
	{
		wspace = _wspace;
	}

	WorkSpace *wspace;
};

#endif

