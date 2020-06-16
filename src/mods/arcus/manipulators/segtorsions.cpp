#include "global.h"
#include "segtorsions.h"
#include "workspace/rotbond.h"
#include "pickers/pickbase.h"

using namespace Maths;

namespace Manipulator
{
	// -----------------------------------------------------------------------------------------------------------
	// Begin: 'SegCoreBBOnlyTC'
	// -----------------------------------------------------------------------------------------------------------

	SegCoreBBOnlyTC::SegCoreBBOnlyTC( SegmentDef& _segdef ) 
		: SegmentDefUser(_segdef),
		m_RotbreakType( _segdef.getBreakType() )
	{
		initialiseRotations();
	}

	void SegCoreBBOnlyTC::initialiseRotations()
	{
		// Set up proxies to workspace to make code access more readable.		
		MoleculeBase& molbase = getWSpace();

		// Clear existing data
		m_RotatePhi_BB.clear();
		m_RotatePsi_BB.clear();
		m_RotateOmega_BB.clear();
		m_CAPos.clear();
		m_BBPos.clear();

		// Initialise PosCache
		int index = -1;
		size_t resIndex = segStartIndex();
		for( size_t i = 0; i < segLength(); i++ )
		{
			index = molbase.res[resIndex].iCA;
			if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iCA");
			m_CAPos.push_back( &molbase.atomxyz(index) );

			index = molbase.res[resIndex].iN;
			if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iN");
			m_BBPos.push_back( &molbase.atomxyz(index) );
			
			// NOTE: The change in the order of the O and C atoms in the container below is IMPORTANT!
			if( m_RotbreakType == SegBreakStart || 
				( m_RotbreakType == SegBreakCentre && resIndex > segCentreIndex() ) 
				)
			{
				// N, O, C
				index = molbase.res[resIndex].iO;
				if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iO");
				m_BBPos.push_back( &molbase.atomxyz(index) );
				
				index = molbase.res[resIndex].iC;
				if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iC");
				m_BBPos.push_back( &molbase.atomxyz(index) );
			}
			else if( m_RotbreakType == SegBreakEnd || 
				( m_RotbreakType == SegBreakCentre && resIndex <= segCentreIndex() ) 
				)
			{
				// N, C, O	
				index = molbase.res[resIndex].iC;
				if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iC");
				m_BBPos.push_back( &molbase.atomxyz(index) );

				index = molbase.res[resIndex].iO;
				if( index == -1 ) THROW(ProcedureException,"Missing atom: Cannot find iO");
				m_BBPos.push_back( &molbase.atomxyz(index) );	
			}
			else
			{
				THROW(NotImplementedException,"Unknown SegBreakType encountered!");
			}

			resIndex++;
		}

		// These are the rotation definitions of the bonds to be rotated for each residues backbone angles
		m_RotatePhi_BB.resize( segLength() );
		m_RotatePsi_BB.resize( segLength() );
		m_RotateOmega_BB.resize( segLength() );

		RotationDefinition_BBOnly *rotDef;

		resIndex = segStartIndex();
		size_t breakIndex = ((segLength()+1) / 2) - 1;
		for( size_t i = 0; i < segLength(); i++ )
		{
			if( m_RotbreakType == SegBreakCentre )
			{
				if( resIndex > segCentreIndex() ) 
				{
					// We are in the second half of the conformer,

					// Do Phi
					rotDef = &m_RotatePhi_BB[i];
					rotDef->conformerIndex = i;
					if( resIndex != 0 )
					{
						rotDef->molBase = &getWSpace();
						rotDef->CAIndexS = breakIndex + 1;
						rotDef->CAIndexE = i - 1; // not including the current CA
						rotDef->BBIndexS = (breakIndex+1) * 3;
						rotDef->BBIndexE = (i * 3) - 1; // Include the pervious residues but not this residues N
						rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
						rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iN );
						rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
						rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iC );
					}

					// Do Psi
					rotDef = &m_RotatePsi_BB[i];
					rotDef->conformerIndex = i;

					// Even if its the last residue, we can still set our "special psi", because we use O not C+1
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = breakIndex + 1;
					rotDef->CAIndexE = i - 1; // not including the current CA
					rotDef->BBIndexS = (breakIndex+1) * 3;
					rotDef->BBIndexE = (i * 3); // Include this residues N
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iC );
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iO );
					// using the 'O' rather than the N of the next residue is non-standard,
					// but simplifying for this program (due to the anchor-2 definition),
					// and valid for idealised geometry systems, i.e. where the peptide group is planar.

					// Do Omega
					rotDef = &m_RotateOmega_BB[i];
					rotDef->conformerIndex = i;
					if( resIndex != 0 )
					{
						rotDef->molBase = &getWSpace();
						rotDef->CAIndexS = breakIndex + 1;
						rotDef->CAIndexE = i - 1; // not including the current CA
						rotDef->BBIndexS = (breakIndex+1) * 3;
						rotDef->BBIndexE = (i * 3) - 2; // the current residue, not including N
						rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex-1].iO );
						rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
						rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iN );
						rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					}
				}
				else
				{
					// Do Phi
					rotDef = &m_RotatePhi_BB[i];
					rotDef->conformerIndex = i;
					if( resIndex != 0 )
					{
						rotDef->molBase = &getWSpace();
						rotDef->CAIndexS = i + 1; // not including the current CA
						rotDef->CAIndexE = breakIndex;
						rotDef->BBIndexS = (i * 3) + 1; // the current residue, not including N
						rotDef->BBIndexE = ((breakIndex+1) * 3) - 1;
						rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
						rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iN );
						rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
						rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iC );
					}

					// Do Psi
					rotDef = &m_RotatePsi_BB[i];
					rotDef->conformerIndex = i;

					// No 'if( resIndex != endRes );' Even if its the last residue, we can still set our 
					// special psi", because we use O not C+1
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = i + 1; // not including the current CA
					rotDef->CAIndexE = breakIndex;
					rotDef->BBIndexS = (i * 3) + 2; // the current residue, not including N or the C
					rotDef->BBIndexE = ((breakIndex+1) * 3) - 1;
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iC );
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iO );
					// using the 'O' rather than the N of the next residue is non-standard,
					// but simplifying for this program (due to the anchor-2 definition),
					// and valid for idealised geometry systems, i.e. where the peptide group is planar.

					// Do Omega
					rotDef = &m_RotateOmega_BB[i];
					rotDef->conformerIndex = i;
					if( resIndex != 0 )
					{
						rotDef->molBase = &getWSpace();
						rotDef->CAIndexS = i; // including the current CA
						rotDef->CAIndexE = breakIndex;
						rotDef->BBIndexS = (i * 3) + 1; // the current residue, not including N
						rotDef->BBIndexE = ((breakIndex+1) * 3) - 1;
						rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex-1].iO );
						rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
						rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iN );
						rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					}
				}	
			}
			else if( m_RotbreakType == SegBreakStart ) 
			{
				// Do Phi
				rotDef = &m_RotatePhi_BB[i];
				rotDef->conformerIndex = i;
				if( resIndex != 0 )
				{
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = 0;
					rotDef->CAIndexE = i - 1; // not including the current CA
					rotDef->BBIndexS = 0;
					rotDef->BBIndexE = (i * 3) - 1; // Include the pervious residues but not this residues N
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iC );
				}

				// Do Psi
				rotDef = &m_RotatePsi_BB[i];
				rotDef->conformerIndex = i;

				// Even if its the last residue, we can still set our "special psi", because we use O not C+1
				rotDef->molBase = &getWSpace();
				rotDef->CAIndexS = 0;
				rotDef->CAIndexE = i - 1; // not including the current CA
				rotDef->BBIndexS = 0;
				rotDef->BBIndexE = (i * 3); // Include this residues N
				rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iN );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iC );
				rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iO );
				// using the 'O' rather than the N of the next residue is non-standard,
				// but simplifying for this program (due to the anchor-2 definition),
				// and valid for idealised geometry systems, i.e. where the peptide group is planar.

				// Do Omega
				rotDef = &m_RotateOmega_BB[i];
				rotDef->conformerIndex = i;
				if( resIndex != 0 )
				{
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = 0;
					rotDef->CAIndexE = i - 1; // not including the current CA
					rotDef->BBIndexS = 0;
					rotDef->BBIndexE = (i * 3) - 2; // the current residue, not including N
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex-1].iO );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
				}
			}
			else if( m_RotbreakType == SegBreakEnd )
			{
				// Do Phi
				rotDef = &m_RotatePhi_BB[i];
				rotDef->conformerIndex = i;
				if( resIndex != 0 )
				{
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = i + 1; // not including the current CA
					rotDef->CAIndexE = segLength() - 1;
					rotDef->BBIndexS = (i * 3) + 1; // the current residue, not including N
					rotDef->BBIndexE = (segLength() * 3) - 1;
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iC );
				}

				// Do Psi
				rotDef = &m_RotatePsi_BB[i];
				rotDef->conformerIndex = i;

				// No 'if( resIndex != endRes );' Even if its the last residue, we can still set our 
				// special psi", because we use O not C+1
				rotDef->molBase = &getWSpace();
				rotDef->CAIndexS = i + 1; // not including the current CA
				rotDef->CAIndexE = segLength() - 1;
				rotDef->BBIndexS = (i * 3) + 2; // the current residue, not including N or the C
				rotDef->BBIndexE = (segLength() * 3) - 1;
				rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex ].iN );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iC );
				rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iO );
				// using the 'O' rather than the N of the next residue is non-standard,
				// but simplifying for this program (due to the anchor-2 definition),
				// and valid for idealised geometry systems, i.e. where the peptide group is planar.

				// Do Omega
				rotDef = &m_RotateOmega_BB[i];
				rotDef->conformerIndex = i;
				if( resIndex != 0 )
				{
					rotDef->molBase = &getWSpace();
					rotDef->CAIndexS = i; // including the current CA
					rotDef->CAIndexE = segLength() - 1;
					rotDef->BBIndexS = (i * 3) + 1; // the current residue, not including N
					rotDef->BBIndexE = (segLength() * 3) - 1;
					rotDef->atom1 = &molbase.atomxyz( molbase.res[resIndex-1].iO );
					rotDef->atom2 = &molbase.atomxyz( molbase.res[resIndex-1].iC );
					rotDef->atom3 = &molbase.atomxyz( molbase.res[resIndex ].iN );
					rotDef->atom4 = &molbase.atomxyz( molbase.res[resIndex ].iCA );
				}
			}
			else
			{
				THROW(NotImplementedException,"Unknown SegBreakType encountered!");
			}
			resIndex++;
		}
	}

	void SegCoreBBOnlyTC::performRotation( RotationDefinition_BBOnly &def, double desiredAngle )
	{
		if( def.CAIndexS < 0 && def.BBIndexS < 0 ) 
			return;

		// Define an axis around which to rotate
		dvector axis;
		axis.setTo( *def.atom3 );
		axis.sub( *def.atom2 );

		// and now ... what rotation matrix is required?
		matrix3x3 rmat;
		rmat.setToAxisRot( axis, desiredAngle - def.getCurrentTorsionAngle() ); 
		
		// Cache the translation locally to avoid pointer lookup in the loop
		dvector trans;
		trans.setTo( *def.atom2 );

		// Loop over CA atoms
		for( int i = def.CAIndexS; i <= def.CAIndexE; i++)
		{			
			Maths::dvector pos( *m_CAPos[i] );
			pos.sub( trans ); // translate to axis origin
			m_CAPos[i]->setTo( 
				( rmat.r[0][0] * pos.x + rmat.r[0][1] * pos.y + rmat.r[0][2] * pos.z ) + trans.x,
				( rmat.r[1][0] * pos.x + rmat.r[1][1] * pos.y + rmat.r[1][2] * pos.z ) + trans.y,
				( rmat.r[2][0] * pos.x + rmat.r[2][1] * pos.y + rmat.r[2][2] * pos.z ) + trans.z
				);
		}
		
		// Loop over BB atoms
		for( int i = def.BBIndexS; i <= def.BBIndexE; i++)
		{			
			Maths::dvector pos = *m_BBPos[i];
			pos.sub( trans ); // translate to axis origin
			m_BBPos[i]->setTo( 
				( rmat.r[0][0] * pos.x + rmat.r[0][1] * pos.y + rmat.r[0][2] * pos.z ) + trans.x,
				( rmat.r[1][0] * pos.x + rmat.r[1][1] * pos.y + rmat.r[1][2] * pos.z ) + trans.y,
				( rmat.r[2][0] * pos.x + rmat.r[2][1] * pos.y + rmat.r[2][2] * pos.z ) + trans.z
				);
		}
	}

	// -----------------------------------------------------------------------------------------------------------
	// End: 'SegCoreBBOnlyTC'
	// -----------------------------------------------------------------------------------------------------------


	// -----------------------------------------------------------------------------------------------------------
	// Begin: 'SegWholeBackboneTC'
	// Used as a base class for seg manipulation classes capable of working on the entrire backbone, performing
	// rotations that manipulate all atoms
	// -----------------------------------------------------------------------------------------------------------

	WholeBackboneTC::WholeBackboneTC() :
		m_InnerMolBase( NULL )
	{
	}

	void WholeBackboneTC::setMolecule( MoleculeBase& _MolBase, SegmentDef& _Range )
	{
		m_InnerMolBase = &_MolBase;
		m_ResRange = _Range;
		initialiseRotations( _Range.getBreakType() );
	}

	void WholeBackboneTC::setMolecule( MoleculeBase& _MolBase, PickResidueRange& _Range, SegBreakType _breakType )
	{
		m_InnerMolBase = &_MolBase;
		m_ResRange = _Range;
		initialiseRotations( _breakType );
	}

	void WholeBackboneTC::initialiseRotations(SegBreakType _breakType)
	{
		ASSERT( m_InnerMolBase != NULL, ProcedureException, "setMolecule() has not been called for WholeBackboneTC");
		
		// get extra initial info from the validated m_ResRange
		MoleculeBase& molbase = *m_InnerMolBase;

		int length = m_ResRange.getNRes();

		// **VERY** important to reinit the whole array - we dont initialise all below
		// the default constructor **HAS** to be called.
		m_RotatePhi_AA.clear();
		m_RotatePsi_AA.clear();
		m_RotateOmega_AA.clear();

		// These are the rotation definitions of the bonds to be rotated for each residues backbone angles
		m_RotatePhi_AA.resize( length );
		m_RotatePsi_AA.resize( length );
		m_RotateOmega_AA.resize( length );

		int residueIndex = -1;
		RotationDefinition_AllAtom *rotDef;

		int outScope = -1;
		int inScope = -1;

		if( _breakType == SegBreakEnd )
		{
			inScope = length; // Disabled
			outScope = length; // the range is the entire loop
		}
		else if( _breakType == SegBreakStart )
		{
			inScope = 0; // the range is the entire loop
			outScope = 0; // Disabled
		}
		else if( _breakType == SegBreakCentre )
		{
			// Half each
			inScope = outScope = ((length+1) / 2);
		}
		else
		{
			THROW(NotImplementedException,"Unknown SegBreakType encountered!");
		}

		int segStartAtomIndex = molbase.res[ m_ResRange.getStartResIndex() + inScope ].ifirst; // the first residues first atom.	
		for( int i = inScope; i < length; i++ )
		{
			residueIndex = m_ResRange.getStartResIndex() + i;

			// Do Phi
			rotDef = &m_RotatePhi_AA[i];
			rotDef->conformerIndex = i;
			if( residueIndex != m_ResRange.getStartResIndex() )
			{
				rotDef->molBase = &molbase;
				rotDef->startAtomIndex = segStartAtomIndex;
				rotDef->endAtomIndex = molbase.res[residueIndex].iN - 1;
				rotDef->arseAtom = molbase.res[residueIndex].iH;
				rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex-1].iC );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex ].iN );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex ].iCA );
				rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex ].iC );
			}

			// Do Psi
			rotDef = &m_RotatePsi_AA[i];
			rotDef->conformerIndex = i;

			// Even if its the last residue, we can still set our "special psi", because we use O not C+1
			// for the torsion definition.
			rotDef->molBase = &molbase;
			rotDef->startAtomIndex = segStartAtomIndex;
			rotDef->endAtomIndex = molbase.res[residueIndex ].iC - 1;
			rotDef->arseAtom = -1;
			rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex].iN );
			rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex].iCA );
			rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex].iC );
			rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex].iO );
			// using the 'O' rather than the N of the next residue is non-standard,
			// but simplifying for this program (due to the anchor-2 definition),
			// and valid for idealised geometry systems, i.e. where the peptide group is planar.

			// DO OMEGA
			rotDef = &m_RotateOmega_AA[i];
			rotDef->conformerIndex = i;
			if( residueIndex != m_ResRange.getStartResIndex() )
			{
				rotDef->molBase = &molbase;
				rotDef->startAtomIndex = segStartAtomIndex;
				rotDef->endAtomIndex = molbase.res[residueIndex-1].iC - 1;
				if( i != 0 ) 
				{
					rotDef->arseAtom = molbase.res[residueIndex-1].iO;
				}
				else
				{
					rotDef->arseAtom = -1;
				}
				rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex-1].iO );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex-1].iC );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex ].iN );
				rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex ].iCA );
			}
		}

		int segEndAtomIndex = molbase.res[ m_ResRange.getStartResIndex() + outScope - 1 ].ilast; // the last residues last atom.	
		for( int i = 0; i < outScope; i++ )
		{
			residueIndex = m_ResRange.getStartResIndex() + i;

			// Do Phi
			rotDef = &m_RotatePhi_AA[i];
			rotDef->conformerIndex = i;
			if( residueIndex != 0 )
			{
				rotDef->molBase = &molbase;
				rotDef->startAtomIndex = molbase.res[residueIndex].iCA + 1; // the current residue. +1 as we are not including N, H or CA
				rotDef->endAtomIndex = segEndAtomIndex;
				rotDef->arseAtom = -1;
				rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex-1].iC );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex ].iN );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex ].iCA );
				rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex ].iC );
			}

			// Do Psi
			rotDef = &m_RotatePsi_AA[i];
			rotDef->conformerIndex = i;
			if( residueIndex != (molbase.nResidues() - 1) )
			{
				rotDef->molBase = &molbase;
				rotDef->startAtomIndex = molbase.res[residueIndex].iO; // the current residue's O is the only atom in the residue that is rotated by psi, all of the next residue must be rotated.
				rotDef->endAtomIndex = segEndAtomIndex;
				rotDef->arseAtom = -1;
				rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex ].iN );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex ].iCA );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex ].iC );
				rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex ].iO );
				// using the 'O' rather than the N of the next residue is non-standard,
				// but simplifying for this program (due to the anchor-2 definition),
				// and valid for idealised geometry systems, i.e. where the peptide group is planar.
			}

			// DO OMEGA
			rotDef = &m_RotateOmega_AA[i];
			rotDef->conformerIndex = i;
			if( residueIndex != 0 )
			{
				rotDef->molBase = &molbase;
				rotDef->startAtomIndex = molbase.res[residueIndex].iN + 1; // not the current residues N atom
				rotDef->endAtomIndex = segEndAtomIndex;
				rotDef->arseAtom = -1;
				rotDef->atom1 = &molbase.atomxyz( molbase.res[residueIndex-1].iO );
				rotDef->atom2 = &molbase.atomxyz( molbase.res[residueIndex-1].iC );
				rotDef->atom3 = &molbase.atomxyz( molbase.res[residueIndex ].iN );
				rotDef->atom4 = &molbase.atomxyz( molbase.res[residueIndex ].iCA );
			}
		}

		//validateAll();
	}

	void WholeBackboneTC::validateAll() const
	{
		// Assert that all the torsions are valid ...
		size_t length = (size_t) m_ResRange.getNRes();
		for( size_t i = 0; i < length; i++ )
		{
			if( !m_RotatePhi_AA[i].isValid() ||
				!m_RotatePsi_AA[i].isValid() ||
				!m_RotateOmega_AA[i].isValid() )
			{
				THROW(ProcedureException,"Invalid angle rotation-group!");
			}
		}
	}

	void WholeBackboneTC::setBackBonePhiPsi_NoValidate( double _Phi, double _Psi )
	{
		ASSERT( m_InnerMolBase != NULL, ProcedureException, "Uninitialised: setMolecule() has not been called for WholeBackboneTC");
		size_t length = (size_t) m_ResRange.getNRes();
		for( size_t i = 0; i < length; i++ )
		{
			m_RotatePhi_AA[i].performRotation(_Phi);
			m_RotatePsi_AA[i].performRotation(_Psi+Maths::MathConst::PI);
		}
	}

	void WholeBackboneTC::setBackBonePhiPsi( double _Phi, double _Psi )
	{
		ASSERT( m_InnerMolBase != NULL, ProcedureException, "Uninitialised: setMolecule() has not been called for WholeBackboneTC");
		size_t length = (size_t) m_ResRange.getNRes();
		for( size_t i = 0; i < length; i++ )
		{
			if( m_RotatePhi_AA[i].isValid() ) 
			{
				m_RotatePhi_AA[i].performRotation(_Phi);
			}
			if( m_RotatePsi_AA[i].isValid() ) 
			{
				m_RotatePsi_AA[i].performRotation(_Psi+Maths::MathConst::PI); // 'PI' due to differing atom definitions to the standard convention
			}
		}
	}

	void WholeBackboneTC::setBackBoneOmega( double _Omg )
	{
		ASSERT( m_InnerMolBase != NULL, ProcedureException, "setMolecule() has not been called for WholeBackboneTC");
		size_t length = (size_t) m_ResRange.getNRes();
		for( size_t i = 0; i < length; i++ )
		{
			if( m_RotateOmega_AA[i].isValid() ) 
			{
				m_RotateOmega_AA[i].performRotation(_Omg+Maths::MathConst::PI); // 'PI' due to differing atom definitions to the standard convention
			}
		}
	}

	//// -----------------------------------------------------------------------------------------------------------
	//// End: 'SegWholeBackboneTC'
	//// -----------------------------------------------------------------------------------------------------------

	// -----------------------------------------------------------------------------------------------------------
	// Begin: 'SegSidechainTC'
	// Used as a base class for seg manipulation classes capable of working on sidechains.
	// There is a protected internal class 'RotBondCache' that is used to store which of the RotBonds we are using.
	// -----------------------------------------------------------------------------------------------------------

	SegSidechainTC::RotBondCache::RotBondCache()
	{
	}

	// calculates the cache for a given sidechain
	void SegSidechainTC::RotBondCache::calcCache( const WorkSpace& _WSpace, int sidechainResNum )
	{
		// Proxies
		int startAtom = _WSpace.res[sidechainResNum].ifirst;
		int endAtom = _WSpace.res[sidechainResNum].ilast;
		const RotBond& rotbond = _WSpace.rotbond();

		// Clear the current cache
		m_RotBondIndexes.clear();

		// Assign the cache ...
		for( size_t i = 0; i < rotbond.size(); i++)
		{
			if( // our bond is a single bond within the current sidechain ...
				(rotbond[i].Type == RotatableBond::Single) &&
				(rotbond[i].i >= startAtom )&&
				(rotbond[i].i <= endAtom )&&
				(rotbond[i].j >= startAtom )&&
				(rotbond[i].j <= endAtom )
				)
			{
				m_RotBondIndexes.push_back(i);
			}
		}
	}

	// calculate which rotatable bonds we are going to use, and add them to a cache.
	void SegSidechainTC::initialiseRotations()
	{
		SegmentDef& segdef = getSegDef();
		int segLength = segdef.getNRes();		

		m_RotBonds.clear();
		m_RotBonds.resize(segLength);		

		for( int i = 0; i < segLength; i++ )
		{
			m_RotBonds[i].calcCache( getWSpace(), segdef.getStartResIndex() + i );
		}
	}

	// -----------------------------------------------------------------------------------------------------------
	// End: 'SegSidechainTC'
	// -----------------------------------------------------------------------------------------------------------

	// -----------------------------------------------------------------------------------------------------------
	// Begin: 'SegTorsionManipulator'
	// -----------------------------------------------------------------------------------------------------------

	SegTorsionManipulator::SegTorsionManipulator( SegmentDef &_segdef ) :
		SegmentDefUser(_segdef),
		SegSidechainTC(_segdef),
		SegWholeBackboneTC(_segdef)
	{
	}

	// ---------------------------------------------------
	// End: 'SegTorsionManipulator'
	// ---------------------------------------------------
}

