#ifndef __SCRACTH_JON_BASICS_H
#define __SCRACTH_JON_BASICS_H

#include "global.h"
#include "forcefields/forcefield.h"
#include "protocols/montecarlo.h"
#include "workspace/workspace.fwd.h"

Physics::Forcefield createffs( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createffts( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createffVac( WorkSpace& wspace, bool useBreakableFF = false, bool summary = true);
Physics::Forcefield createff(WorkSpace& wspace, bool useBreakableFF = false, double dielec = 1.0, bool summary = true);

double getMeRMS( const std::vector<Maths::dvector>& native, const std::vector<Maths::dvector>& conformer );

#endif

