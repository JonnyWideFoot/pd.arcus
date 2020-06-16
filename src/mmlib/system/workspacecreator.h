#ifndef __WORKSPACE_CREATOR_H
#define __WORKSPACE_CREATOR_H

#include "workspace/workspace.fwd.h"
namespace Sequence
{
	class PD_API BioSequence;
}

class PD_API System;
class PD_API FFParamSet;






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
class PD_API WorkspaceCreatorBase
{
public:
	friend class PD_API WorkSpace;
	WorkspaceCreatorBase(System &_system);

	// Public accessors
	const FFParamSet &ffps() const;

protected:
	Sequence::BioSequence& WSCall_GetSeq(WorkSpace &wspace) const;
	void WSCall_Allocate(WorkSpace &wspace, int _natom) const;
	void WSCall_ReInitialiseAll(WorkSpace &wspace) const;
	void WSCall_Append(WorkSpace &wspace, const System& sysspec) const;

	void create(WorkSpace &wspace) const;
	virtual void allocateParticles(WorkSpace &wspace) const = 0;

	void setAtomPositionRedirection(WorkSpace &wspace) const;

	System *system;
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
class PD_API WSfull: public WorkspaceCreatorBase
{
public:
	friend class PD_API WorkSpace;
	WSfull(System &system):WorkspaceCreatorBase(system) {}

protected:
	virtual void allocateParticles(WorkSpace &wspace) const;
};

#endif

