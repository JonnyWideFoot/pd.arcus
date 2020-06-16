#ifndef __SEGBUILDER_H
#define __SEGBUILDER_H

#include "system/rebuilder.h"

class MoleculeBase;

// ---------------------
//  Rebuilder Functions
// ---------------------

bool rebuildExtdFromStart(MoleculeBase& mol, int start, int end);
bool rebuildExtdFromEnd(MoleculeBase& mol, int start, int end);

/// Returns a vector of ResidueRanges that require rebuild
/// (i.e. whole sections of residues where isRebuildRequired() == true)
std::vector<PickResidueRange> findMissingSegments(MoleculeBase& mol);

// -------------------
//  Rebuilder Classes
// -------------------

//-------------------------------------------------
//
/// \brief
/// Builds missing whole sections from their closest anchor residue on their polymer.
/// Then rebuilds remaining missing atoms using forcefield definitions.
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author  Jon Rea
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
class SegmentRebuilder : public RebuilderBase
{
public:
	SegmentRebuilder();
	SegmentRebuilder(bool _DoDualBranch);
	virtual bool invokeBuild(MoleculeBase& mol, Verbosity::Type verbose);
	void setDualBranch( bool _DualBranchMode ) { m_DualBranch = _DualBranchMode; }

	// Public flags

	/// Defines that once the sections have been extended and stitched, they will once again be
	/// marked as 'isRebuildRequired'. This is the default, allowing more advanced
	/// workspace/protocol-based procedures to be used later on with full forcefields.
	bool DownstreamFullSectionRebuild;

	const std::vector<PickResidueRange>& getPreviousBuildSections() const { return m_Seg; }

protected:
	std::vector<PickResidueRange> m_Seg;
	bool buildMissingSegments(MoleculeBase& mol);
	bool buildFromStart(MoleculeBase& mol, int start, int end);
	bool buildFromEnd(MoleculeBase& mol, int start, int end);
	void unsetSectionBuildFlag(MoleculeBase& mol);
	bool m_DualBranch;
};

#endif

