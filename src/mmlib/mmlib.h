#ifndef __MMLIB_H
#define __MMLIB_H

// Should contain all the mmlib headers important for modules, python wrapping etc..

#include "mmlib/global.h"
#include "mmlib/interface.h"
#include "mmlib/object.h"
#include "mmlib/primitives.h"
#include "mmlib/exception.h"

#include "mmlib/maths/maths.h"
#include "mmlib/maths/maths_matrix3x3.h"
#include "mmlib/maths/maths_vector.h"
#include "mmlib/maths/graphtheory.h"
#include "mmlib/maths/fastrandom.h"

#include "mmlib/maths/histogram.h"
#include "mmlib/tools/statclock.h"
#include "mmlib/tools/streamwriter.h"
#include "mmlib/tools/stringbuilder.h"

#include "mmlib/system/fundamentals.h"
#include "mmlib/system/molecule.h"
#include "mmlib/system/system.h"
#include "mmlib/system/genpolymer.h"
#include "mmlib/pickers/pickbase.h"
#include "mmlib/pickers/basicpickers.h"
#include "mmlib/pickers/pickfromfile.h"
#include "mmlib/system/rebuilder.h"

#include "mmlib/library/backbonetorsions.h"
#include "mmlib/library/angleset.h"
#include "mmlib/library/mapper.h"
#include "mmlib/library/nameconv.h"
#include "mmlib/library/rotamerlib.h"
#include "mmlib/library/rotamer_dunbrack.h"
#include "mmlib/library/rotamer_shetty.h"

#include "mmlib/sequence/sequence.h"
#include "mmlib/sequence/alignment.h"

#include "mmlib/fileio/outtra.h"
#include "mmlib/fileio/intra.h"
#include "mmlib/fileio/infile.h"
#include "mmlib/fileio/pdb.h"
#include "mmlib/fileio/tra.h"
#include "mmlib/fileio/trablocks.h"
#include "mmlib/fileio/tratypes.h"
#include "mmlib/fileio/psfdcd.h"

#include "mmlib/workspace/componentbase.h"
#include "mmlib/workspace/space.h"
#include "mmlib/workspace/bondorder.h"
#include "mmlib/workspace/neighbourlist.h"
#include "mmlib/workspace/rotbond.h"
#include "mmlib/workspace/workspace.h"
#include "mmlib/workspace/snapshot.h"
#include "mmlib/workspace/pospointer.h"

#include "mmlib/sequence/sequence.h"
#include "mmlib/sequence/alignment.h"

#include "mmlib/library/residues.h"
#include "mmlib/library/nameconv.h"

#include "mmlib/forcefields/ffparam.h"
#include "mmlib/forcefields/forcefield.h"
#include "mmlib/forcefields/ffbonded.h"
#include "mmlib/forcefields/gbff.h"
#include "mmlib/forcefields/lcpo.h"

#include "mmlib/forcefields/restraintbase.h"
#include "mmlib/forcefields/restraint_positional.h"
#include "mmlib/forcefields/restraint_internal.h"
#include "mmlib/forcefields/restraint_torsional.h"
#include "mmlib/forcefields/restraint_atomdist.h"
#include "mmlib/forcefields/restraint_native_contact.h"

#include "mmlib/forcefields/ffcustom.h"
#include "mmlib/forcefields/ffsoftvdw.h"
#include "mmlib/forcefields/nonbonded.h"
#include "mmlib/forcefields/nonbonded_ti.h"
#include "mmlib/forcefields/nonbonded_ti_linear.h"
#include "mmlib/forcefields/breakablebonded.h"

#include "mmlib/manipulators/movebase.h"
#include "mmlib/manipulators/rotamer_applicatorbase.h"
#include "mmlib/manipulators/basicmoves.h"
#include "mmlib/manipulators/rotamer_scwrl.h"

#include "mmlib/filters/filterbase.h"
#include "mmlib/filters/basicfilters.h"

#include "mmlib/protocols/protocolbase.h"
#include "mmlib/protocols/energy.h"
#include "mmlib/protocols/minimise.h"
#include "mmlib/protocols/temperature.h"
#include "mmlib/protocols/md.h"
#include "mmlib/protocols/rerun.h"
#include "mmlib/protocols/remd.h"
#include "mmlib/protocols/torsionalminimisation.h"
#include "mmlib/protocols/scpack.h"
#include "mmlib/protocols/dualffminimiser.h"
#include "mmlib/protocols/montecarlo.h"

#include "mmlib/monitors/basicmonitors.h"
#include "mmlib/monitors/umbrella.h"
#include "mmlib/monitors/monexp.h"
#include "mmlib/protocols/nmode.h"

#endif


