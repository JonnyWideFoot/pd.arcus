#ifndef __SC_PACK
#define __SC_PACK

#include <string>

#include "workspace/workspace.fwd.h"
#include "forcefields/forcefield.fwd.h"

namespace Library
{
	class RotamerLibrary;
}

/// Defines a monte carlo with minimiation to pack sidechains only, unaffecting the backbone
void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const Library::RotamerLibrary& rotLib );
void PD_API MCPackSideChains( WorkSpace& wspace, Physics::Forcefield& ffs, Physics::Forcefield& ff, const std::string& rotLibPath );

#endif

