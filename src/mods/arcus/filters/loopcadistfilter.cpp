#include "global.h"

#include "pickers/basicpickers.h"
#include "workspace/workspace.h"
#include "system/proximitygrid.h"

#include "loopcadistfilter.h"

LoopCADistFilter::LoopCADistFilter() 
	: SqrCAMinReach(6.32*6.32), // See PLOP method paper for why 6.32
	m_NOffset(0),
	m_COffset(0)
{
}

void LoopCADistFilter::setTo( const ProximityGrid& _StaticCAGrid, const MoleculeBase& mol, const PickResidueRange& picker, size_t _NOffset, size_t _COffset )
{
	m_Mol = &mol;
	m_StaticCAGrid = &_StaticCAGrid;
	m_Picker = picker;
	m_NOffset = _NOffset;
	m_COffset = _COffset;
}

bool LoopCADistFilter::passes()
{
	size_t from = m_Picker->getStartResIndex() + m_NOffset;
	size_t to = m_Picker->getEndResIndex() - m_COffset;
	for( size_t i = from; i <= to; i++ )
	{
		bool foundOne = false;
		const Residue& res = m_Mol->res[i];

		size_t to_j = res.ilast;
		for( size_t j = res.ifirst; j <= to_j; j++ )
		{
			const GridPoint* list = NULL;
			const Particle& particleJ = m_Mol->atom[j];
			if( !particleJ.isCoreBackbone() )
				continue; // they just dont count, baby :-D

			const Maths::dvector& atomJPos = particleJ.pos();
			m_StaticCAGrid->atomList( atomJPos, list );

			if( list == NULL || list->size() == 0 )
			{
				// There is NOTHING (that is static) even remotely close to this atom in space.
				return false;
			}

			foundOne = false;
			for( size_t k = 0; k < list->size(); k++ )
			{
				int gridAtom = list->atomIndexes[k];
				const Maths::dvector& atomKPos = m_Mol->atom[gridAtom].pos();
				double sqrDist = atomJPos.sqrdist( atomKPos );
				if( sqrDist < SqrCAMinReach )
				{
					// Cool, we are close enough to the protein body for this residue :-D
					// Continue to check the remaining residues...
					foundOne = true;
					break;
				}
			}
			if( !foundOne )
				return false;			
		}
	}
	// All loop residue backbones have a contact < 6.32A from a heavy atom in the protein body
	return true;
}

std::string LoopCADistFilter::reason()
{
	StringBuilder sb;
	size_t from = m_Picker->getStartResIndex() + 1;
	size_t to = m_Picker->getEndResIndex() - 1;
	for( size_t i = from; i <= to; i++ )
	{
		bool foundOne = false;
		const Residue& res = m_Mol->res[i];

		size_t to_j = res.ilast;
		for( size_t j = res.ifirst; j <= to_j; j++ )
		{
			const GridPoint* list = NULL;
			const Particle& particleJ = m_Mol->atom[j];
			if( !particleJ.isCoreBackbone() )
				continue; // they just dont count, baby :-D

			const Maths::dvector& atomJPos = particleJ.pos();
			m_StaticCAGrid->atomList( atomJPos, list );

			if( list == NULL || list->size() == 0 )
			{
				// There is NOTHING (that is static) even remotely close to this atom in space.
				if( sb.size() == 0 ) sb.setFormat("Residue %d, no backbone atoms have a contact <6.32A")(i);
				return sb.toString();
			}

			foundOne = false;
			for( size_t k = 0; k < list->size(); k++ )
			{
				int gridAtom = list->atomIndexes[k];
				const Maths::dvector& atomKPos = m_Mol->atom[gridAtom].pos();
				double sqrDist = atomJPos.sqrdist( atomKPos );
				if( sqrDist < SqrCAMinReach )
				{
					// Cool, we are close enough to the protein body for this residue :-D
					// Continue to check the remaining residues...
					foundOne = true;
					break;
				}
			}
			if( !foundOne )
			{
				sb.setFormat("Residue %d, backbone atom %d, %s has no contact <6.32A")(i)(j)(particleJ.pdbname.c_str());			
				return sb.toString();
			}
		}
	}
	// All loop residue backbones have a contact < 6.32A from a heavy atom in the protein body
	return "All loop backbone atoms have protein body contact < 6.32A";	
}

