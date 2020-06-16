#ifndef __BREAKABLE_BONDED_H
#define __BREAKABLE_BONDED_H

#include "tools/vector.h"
#include "primitives.h"
#include "ffbonded.h"
#include "ffcustom.h"

namespace Physics
{
	//-------------------------------------------------
	/// \brief  Abstract Base Class or 'Interface' to a forcefield that support broken bonds.
	/// \author Jon Rea 
	class PD_API IBondBreakable
	{
	public:
		virtual void createBreak( int _AtomIndexA, int _AtomIndexB ) = 0;
		virtual void clearBreaks() = 0;
	};


	//-------------------------------------------------
	/// \brief  A bonded-component forcefield that allows the removal of selected bonded terms.
	/// \details For a given broken bond, the bond angle torsion and improper terms are removed from
	/// the calculation. This is obviously only viable in certain circumstances.
	/// \author Jon Rea 
	class PD_API FF_BreakableBonded : public FF_Bonded, public IBondBreakable
	{
		struct BreakDef : public IndexPair
		{
			BreakDef() 
			{ 
				valid = false; 
				i = -1;
				j = -1;
			}
			BreakDef( int _i, int _j ) 
			{ 
				valid = false; 
				i = _i;
				j = _j;
			}
			bool valid; ///< Represents an actual chemical bond
			std::vector<int> atomList; ///< stores the indexes from one side of this break definition
		};

	public:
		FF_BreakableBonded(WorkSpace &newwspace);

		virtual FF_BreakableBonded* clone() const { return new FF_BreakableBonded(*this); }

		virtual void createBreak( int _AtomIndexA, int _AtomIndexB );
		virtual void clearBreaks();

	protected:

		void validateBonds();
		virtual void setup();

		// Hold the break definitions
		std::vector<BreakDef> m_BrokenBonds;

		bool shouldRemove( const Physics::Bond bond ) const;
		bool shouldRemove( const Physics::Angle angle ) const;
		bool shouldRemove( const Physics::Torsion torsion ) const;
	};
}

#endif

