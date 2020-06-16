#ifndef __ARCUS_H
#define __ARCUS_H

#include "arcus_refiner.h"

namespace Library
{
	class RotamerLibrary;
}

namespace Protocol
{
	class PD_API Arcus : public ArcusBase, public IO::BTF_Block_Comment_User
	{
	public:
		Arcus( 
			Physics::Forcefield& _ffs, 
			Physics::Forcefield& _ff, 
			const Library::AngleSet& _as, 
			const Library::RotamerLibrary& _rotLib );

		virtual Arcus* clone() const 
		{ 
			return new Arcus(*this); 
		}

		virtual int runcore();

		size_t ProduceNJoinedBranches;

		// Assert that basic filters are not removing the native!
		void AssertCA6_32Filter(); ///< Probably a good plan to assert this works for each native of interest...
		void AssertSegDistFilter(); ///< Probably a good plan to assert this works for each native of interest...

		void setDefaultRefiner();
		void setRefiner( ArcusRefineBase& _refiner );

		std::string nameMe; ///< TEMPORARY HACK - Delete me after development!!

	protected:

		virtual int initialise();

		virtual void configLevel1Filters( size_t i );

		void calibrateBuilders();
		void generateValidConformers(); // Stage-1
		void rejoinConformers(); // Stage-2
		void refineConformers(); // Stage-3

		std::vector< LoopSet > m_BuiltSections;

		const Library::RotamerLibrary* m_RotLib;

	private:
		PickAtomRanges m_DynamicAtomsPicker;
		ProximityGrid m_StaticGrid4A; ///< A proximity grid containing all atoms which are not involved in the m_Regions array - i.e. those that dont move. 4.0 angstrom cutoff designed for clash-detection
		ProximityGrid m_StaticGrid6_5A;///< A proximity grid containing all atoms which are not involved in the m_Regions array - i.e. those that dont move. 7.0 angstrom cutoff designed for loop:rigid-body contact detection.

		ArcusRefine_CGMin m_DefaultStage3;
		ArcusRefineBase* m_Stage3Refiner;
	};
}

#endif

