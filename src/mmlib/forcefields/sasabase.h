#ifndef __SASA_BASE_H
#define __SASA_BASE_H

#include <vector>
#include "tools/cloneholder.h" // Provides class member
#include "forcefields/forcefield.h" // provides base class

namespace Physics
{





//-------------------------------------------------
//
/// \brief  
/// The common properties for all surface area forcefields
/// Additional arrays will be required in the parent class for additional per atom properties.
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
	class SurfaceEntry
	{
	public:
		int parentAtomIndex;
		Maths::dvector* parentAtomPos;
		double radius;
		double sigma;
	};






//-------------------------------------------------
//
/// \brief  The base class for any forcefields that deal with surface area
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
	class PD_API FF_SurfaceArea : public ForcefieldBase
	{
	public:
		FF_SurfaceArea(WorkSpace &newwspace);
		virtual ~FF_SurfaceArea();
		virtual FF_SurfaceArea* clone() const = 0;

		virtual void setToDefaults();
		void setForPicker( const PickBase& _ForPicker );
		void setAgainstPicker( const PickBase& _AgainstPicker );

		// Public accessors - obtain the specifics of the previous calculation
		//double sasaFraction( const PickBase& _Range ) const;
		//double SASA( const PickBase& _Range ) const;
		//double sasaFraction( const PickAtomRange& _Range ) const;
		//double SASA( const PickAtomRange& _Range ) const;
		//double atomSASA( int ia ) const;
		//double atomFraction( int ia ) const;
		//double resSASA( int ir ) const;
		//double resFraction( int ir ) const;

	protected:
		CloneHolder<PickBase> m_ForPicker; ///< Dictates which atoms we will calculate a SASA for.
		CloneHolder<PickBase> m_AgainstPicker; ///< Dictates which atoms we will compare the 'for' atoms against for solvent occlusion

		WorkSpace* m_WSpace; ///< The current parent workspace	
		std::vector<size_t> m_AtomIndexes; ///< An array containin the index in m_Atoms for a given atom, or SIZE_T_FAIL if it is not in the list
		std::vector<SurfaceEntry> m_Atoms; ///< Array containing a subset of atoms which point to their parents and contain the relevent parameters
	};
}

#endif

