#include "global.h"
#include "tools/vector.h"
#include "system/molecule.h"
#include "system/genpolymer.h"
#include "manipulators/segtorsions.h"
#include "segbuilder.h"

using namespace Maths;

std::vector<PickResidueRange> findMissingSegments(MoleculeBase& mol)
{
	std::vector<PickResidueRange> segList;

	bool extending = false;
	size_t size = mol.res.size();
	size_t start = SIZE_T_FAIL;

	for( size_t i = 0; i < size; i++ )
	{
		int known = 0;
		for( int j = mol.res[i].ifirst; j <= mol.res[i].ilast; j++ )
		{
			if( !mol.atom[j].isRebuildRequired() )
			{
				known++;
			}
		}
		int atomCount =  mol.res[i].ilast - mol.res[i].ifirst + 1;
		if( atomCount > 3 && known < 3 )
		{
			// Essentially the whole residue isn't known 
			// (you cant even rebuild with <3 atoms)
			if( !extending )
			{
				start = i;
			}
			extending = true;
		}
		else
		{
			// residue is known, we are fine here :-D
			if( extending )
			{
				PickResidueRange segDef(mol,start,i-1);
				segList.push_back(segDef);
			}
			extending = false;
		}
	}

	if( extending )
	{
		PickResidueRange segDef(mol,start,size-1);
		segList.push_back(segDef);
	}

	return segList;
}

bool rebuildExtdFromStart(MoleculeBase& mol, int start, int end)
{
	if( start == 0 )
	{
		printf("ERROR: Cannot 'rebuildExtdFromStart' if the first residue is the starting terminus!");
		return false;
	}

	// Local reference to the molecules atom container
	ParticleStore & atom = mol.atom;
	const ResidueStore& res = mol.res;

	printf("Proceeding to build %d to %d from the initiating anchor\n",start,end);

	// --------------------------------------
	//  STAGE 1: Find required stiching atoms
	// --------------------------------------

	// Find out if subsequent stiching is possible ...
	int superimpositionStart = start-1;
	if( !res[start].param->backlink )
	{
		printf("The stitching residue has no backwards link. Reconnection is impossible!");
		return false;
	}
	int iBackAtom = mol.findParticle(superimpositionStart,res[start].param->backname.c_str());
	if(iBackAtom < 0) {
		printf("ERROR: Cannot find Backward Link Atom in residue %d\n", superimpositionStart);
		return false;
	}
	if( !res[superimpositionStart].param->frwdlink )
	{
		printf("The stitching residue has no forward link. Reconnection is impossible!");
		return false;
	}
	int iFrwdAtom = mol.findParticle(start,res[superimpositionStart].param->frwdname.c_str());
	if(iFrwdAtom < 0) {
		printf("ERROR: Cannot find Forward Link Atom in residue %d\n", start);
		return false;
	}

	// Ensure that the next residue has all the relevent stitching atoms
	printf(" Ensuring anchor residue %d is completely rebuilt.\n",superimpositionStart);
	if( !rebuildMissingAtoms( mol, superimpositionStart, Verbosity::Normal ) )
	{
		printf("Failed anchor residue rebuild!!\n");
		return false;
	}

	// *IMPORTANT* - we must make sure that the forward link is not included in the superimposition list generated below
	// Marking the atom as requiring rebuild will ensure this is the case.
	atom[iFrwdAtom].setRebuildRequired(true);

	std::vector<int> r1;
	std::vector<int> r2;

	// Our first atom is the backward link...
	r1.push_back( iBackAtom );

	int b0, b1, b2 = -1;

	int found_R2 = findLevelTwo( atom, r1, r2 );
	if( found_R2 == 0 )
	{
		// We are doomed; There are no level 2 atoms.
	}
	else if( found_R2 >= 2 )
	{
		b0 = r1[0];
		b1 = r2[0];
		b2 = r2[1];
	}
	else if(findLevelThree( atom, r1, r2, b2 ))
	{
		b0 = r1[0];
		b1 = r2[0];
	}

	// Restore this status - see comment above.
	atom[iFrwdAtom].setRebuildRequired(false);

	// Did we get the required atoms??
	if( b0 == -1 )
	{
		printf("Appropriate stitching atoms can't be found!\n");
		return false;
	}
	else
	{
		printf(" Stitching using atoms:\n");
		atom[b0].info();
		atom[b1].info();
		atom[b2].info();
	}

	// --------------------------------------
	//  STAGE 2: Rebuild Operations
	// --------------------------------------

	// Firstly record the current atom positions of residue 'start-1' - this will be used for
	// superimposition.
	int backupStart = res[superimpositionStart].ifirst;
	int backupEnd = res[superimpositionStart].ilast;
	int backupCount =  backupEnd - backupStart + 1;
	std::vector<dvector> origPos;
	origPos.resize(backupCount);
	for( int i = 0; i < backupCount; i++ )
	{
		origPos[i].setTo(mol.atomxyz(i+backupStart));
	}

	// Polymerise our new atomic structure from the forward and back link definititions
	if( !rebuildAndPolymerise( mol, superimpositionStart, end ) ) return false;

	// --------------------------------------
	//  STAGE 3: Perform stiching operation
	// --------------------------------------

	matrix3x3 rmat;

	superimpose(
		origPos[b0-backupStart],
		origPos[b1-backupStart],
		origPos[b2-backupStart],
		mol.atomxyz(b0),
		mol.atomxyz(b1),
		mol.atomxyz(b2),
		rmat);

	int fromAtom = res[start].ifirst;
	int toAtom = res[end].ilast;
	for( int i = fromAtom; i <= toAtom; i++ )
	{
		mol.atomxyz(i).sub(mol.atomxyz(b0));
		mol.atomxyz(i).mulmat(rmat);
		mol.atomxyz(i).add(origPos[b0-backupStart]);
	}

	// Restore the anchor atoms to their original positions
	for( int i = 0; i < backupCount; i++ )
	{
		mol.atomxyz(i+backupStart).setTo(origPos[i]);
	}

	// Attempt torsional rotation (we assume its a polypeptide at the mo)
	// This call will silently fail if not, but the bond and angle geomatries 
	// will be fine!
	Manipulator::WholeBackboneTC tc;
	PickResidueRange def( mol, start, end );
	tc.setMolecule(mol,def,SegBreakEnd);
	tc.setBackBonePhiPsi( DegToRad(-57.0), DegToRad(-47.0) );
	tc.setBackBoneOmega( DegToRad(180.0) );

	return true;
}

bool rebuildExtdFromEnd(MoleculeBase& mol, int start, int end)
{
	if( end == mol.res.size()-1 )
	{
		printf("ERROR: Cannot 'rebuildExtdFromEnd' if the last residue is the polymer terminus!");
		return false;
	}

	// Local reference to the molecules atom container
	ParticleStore & atom = mol.atom;
	const ResidueStore& res = mol.res;

	printf("Proceeding to build %d to %d from the terminating anchor\n",start,end);

	// --------------------------------------
	//  STAGE 1: Find required stiching atoms
	// --------------------------------------

	// Find out if subsequent stiching is possible ...
	int superimpositionEnd = end+1;
	if( !res[superimpositionEnd].param->backlink )
	{
		printf("The stitching residue has no backward link. Reconnection is impossible!");
		return false;
	}
	int iBackAtom = mol.findParticle(end,res[superimpositionEnd].param->backname.c_str());
	if(iBackAtom < 0) {
		printf("ERROR: Cannot find Backward Link Atom in residue %d\n", end);
		return false;
	}
	if( !res[end].param->frwdlink )
	{
		printf("The stitching residue has no forward link. Reconnection is impossible!");
		return false;
	}
	int iFrwdAtom = mol.findParticle(superimpositionEnd,res[end].param->frwdname.c_str());
	if(iFrwdAtom < 0) {
		printf("ERROR: Cannot find Forward Link Atom in residue %d\n", superimpositionEnd);
		return false;
	}

	// Ensure that the next residue has all the relevent stitching atoms
	printf(" Ensuring anchor residue %d is completely rebuilt.\n",superimpositionEnd);
	if( !rebuildMissingAtoms( mol, superimpositionEnd, Verbosity::Normal ) )
	{
		printf("Failed anchor residue rebuild!!\n");
		return false;
	}

	// *IMPORTANT* - we must make sure that the forward link is not included in the superimposition list generated below
	// Marking the atom as requiring rebuild will ensure this is the case.
	atom[iBackAtom].setRebuildRequired(true);

	std::vector<int> r1;
	std::vector<int> r2;

	// Our first atom is the backward link...
	r1.push_back( iFrwdAtom );

	int b0, b1, b2 = -1;

	int found_R2 = findLevelTwo( atom, r1, r2 );
	if( found_R2 == 0 )
	{
		// We are doomed; There are no level 2 atoms.
	}
	else if( found_R2 >= 2 )
	{
		b0 = r1[0];
		b1 = r2[0];
		b2 = r2[1];
	}
	else if(findLevelThree( atom, r1, r2, b2 ))
	{
		b0 = r1[0];
		b1 = r2[0];
	}

	// Restore this status - see comment above.
	atom[iBackAtom].setRebuildRequired(false);

	// Did we get the required atoms??
	if( b0 == -1 )
	{
		printf("Appropriate stitching atoms can't be found!");
		return false;
	}
	else
	{
		printf(" Stitching using atoms:\n");
		atom[b0].info();
		atom[b1].info();
		atom[b2].info();
	}

	// --------------------------------------
	//  STAGE 2: Rebuild Operations
	// --------------------------------------

	// Firstly record the current atom positions of residue 'start-1' - this will be used for
	// superimposition.
	int backupStart = res[superimpositionEnd].ifirst;
	int backupEnd = res[superimpositionEnd].ilast;
	int backupCount =  backupEnd - backupStart + 1;
	std::vector<dvector> origPos;
	origPos.resize(backupCount);
	for( int i = 0; i < backupCount; i++ )
	{
		origPos[i].setTo(mol.atomxyz(i+backupStart));
	}

	// Polymerise our new atomic structure from the forward and back link definititions
	if( !rebuildAndPolymerise( mol, start, superimpositionEnd ) ) return false;

	// --------------------------------------
	//  STAGE 3: Perform stiching operation
	// --------------------------------------

	matrix3x3 rmat;

	superimpose(
		origPos[b0-backupStart],
		origPos[b1-backupStart],
		origPos[b2-backupStart],
		mol.atomxyz(b0),
		mol.atomxyz(b1),
		mol.atomxyz(b2),
		rmat);

	int fromAtom = res[start].ifirst;
	int toAtom = res[end].ilast;
	for( int i = fromAtom; i <= toAtom; i++ )
	{
		mol.atomxyz(i).sub(mol.atomxyz(b0));
		mol.atomxyz(i).mulmat(rmat);
		mol.atomxyz(i).add(origPos[b0-backupStart]);
	}

	// Restore the anchor atoms to their original positions
	for( int i = 0; i < backupCount; i++ )
	{
		mol.atomxyz(i+backupStart).setTo(origPos[i]);
	}

	Manipulator::WholeBackboneTC tc;

	PickResidueRange def( mol, start, end );
	tc.setMolecule(mol,def,SegBreakStart);
	tc.setBackBonePhiPsi( DegToRad(-57.0), DegToRad(-47.0) );

	// need a second defintion to ensure that the omega angle to the stitching residue is set!
	PickResidueRange def2( mol, start, end + 1 ); 
	tc.setMolecule(mol,def2,SegBreakStart);
	tc.setBackBoneOmega( DegToRad(180.0) );

	return true;
}

SegmentRebuilder::SegmentRebuilder() :
	m_DualBranch(false),
	DownstreamFullSectionRebuild(true)
{
}

SegmentRebuilder::SegmentRebuilder(bool _DoDualBranch) : 
	m_DualBranch(_DoDualBranch),
	DownstreamFullSectionRebuild(true)
{
}

bool SegmentRebuilder::invokeBuild(MoleculeBase& mol, Verbosity::Type verbose)
{
	m_Seg = findMissingSegments(mol);
	if(!buildMissingSegments(mol)) 
	{
		return false;
	}
	if( !rebuildMissingAtoms(mol,verbose) )
	{
		return false;
	}	
	if( DownstreamFullSectionRebuild )
	{
		unsetSectionBuildFlag(mol);
	}	
	return true;
}

void SegmentRebuilder::unsetSectionBuildFlag(MoleculeBase& mol)
{
	for( size_t i = 0; i < m_Seg.size(); i++ )
	{
		// This flag determins that something cleaver, probably a protocol, 
		// down the code execution path will do the full rebuild! 
		// Probably will be at the workspace level and require a 
		// forcefield of some description :-D

		ParticleStore & atom = mol.atom;
		const ResidueStore& res = mol.res;
		for( size_t i = 0; i < m_Seg.size(); i++ )
		{
			size_t start = m_Seg[i].getStartResIndex();
			size_t end = m_Seg[i].getEndResIndex();
			int toAtom = res[end].ilast;
			for( int i = res[start].ifirst; i <= toAtom; i++ )
			{
				atom[i].setRebuildRequired(true);
			}
		}
	}
}

bool SegmentRebuilder::buildMissingSegments(MoleculeBase& mol)
{
	if( m_Seg.size() == 0 ) 
	{
		return true; // nothing to do!
	}

	for( size_t i = 0; i < m_Seg.size(); i++ )
	{
		size_t start = m_Seg[i].getStartResIndex();
		size_t end = m_Seg[i].getEndResIndex();

		if( start == 0 )
		{
			if( !rebuildExtdFromEnd(mol,start,end) ) 
				return false; // the other end doesn't exist, so there is nothing to build on...
		}
		else if ( end == mol.res.size()-1 )
		{
			if( !rebuildExtdFromStart(mol,start,end) ) 
				return false; // the other end doesn't exist, so there is nothing to build on...
		}
		else
		{
			int length = end - start + 1;
			if( length < 2 )
			{
				// Uber short! ... Default to 'FromStart' build Mode
				if( !rebuildExtdFromStart(mol,start,end) ) 
					return false;
			}
			else if( m_DualBranch )
			{
				int midPoint = start + (length / 2);
				int remainder = length % 2;
				if( !rebuildExtdFromStart(mol,start,midPoint + remainder) ) 
					return false;
				if( !rebuildExtdFromEnd(mol,midPoint+1,end)) 
					return false;				
			}
			else 
			{
				// Default to Start build Mode
				if( !rebuildExtdFromStart(mol,start,end) ) 
					return false;
			}
		}
	}

	return true;
}

