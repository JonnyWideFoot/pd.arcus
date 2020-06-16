#ifndef __SEG_TORSIONS_H
#define __SEG_TORSIONS_H

#include <vector>

#include "workspace/rotbond.h" // Required for base class
#include "workspace/segdef.h" // Required for base class
#include "pickers/pickbase.h"

namespace Manipulator
{
	//-------------------------------------------------
	/// \brief  Intended as a base class for ConformerBuilder,
	/// \details
	/// to apply backbone torsions to the particle system within
	/// the seg structure only and only the 4 main backbone atoms ...
	/// Used as a base class for seg manipulation classes capable of working (temporarily) on ONLY the main backbone
	/// atoms, totally disregarding the other atoms. This is obviously only applicable in a limited set of
	/// circumstances.
	/// \author  Jon Rea 
	class PD_API SegCoreBBOnlyTC : public virtual SegmentDefUser
	{
	protected:
		//-------------------------------------------------
		/// \brief  RotationDefinition_BBOnly: Intended for use within SegCoreBBOnlyTC
		/// \details Contains the defintions of which atoms need moving per rotation
		/// \author Jon Rea 
		class RotationDefinition_BBOnly
		{
		public:
			MoleculeBase *molBase; // internal pointer set at the same time that the

			Maths::dvector *atom1; // anchor atom
			Maths::dvector *atom2; // atom of the torsion bond
			Maths::dvector *atom3; // atom of the torsion bond
			Maths::dvector *atom4; // anchor atom

			int conformerIndex;

			// used to rotate just the CA and BB atom caches - the index is dependent on the _RotBackward flag
			// in the parent class. This index marks the midpoint, dependent on wether we are going from the start 
			// or the end of the segment.
			int CAIndexS; ///< -1 to set this class to null status
			int BBIndexS;
			int CAIndexE; ///< -1 to set this class to null status
			int BBIndexE;

			RotationDefinition_BBOnly():
			molBase(NULL),
				conformerIndex(-1),
				CAIndexS(-1),
				BBIndexS(-1),
				CAIndexE(-1),
				BBIndexE(-1),
				atom1(NULL),
				atom2(NULL),
				atom3(NULL),
				atom4(NULL)
			{
			}

			inline double getCurrentTorsionAngle(){ return calcTorsionAngle( *atom1, *atom2, *atom3, *atom4 ); }
		};

	public:
		SegCoreBBOnlyTC( SegmentDef& _segdef );

	protected:
		void initialiseRotations();

		void performRotation( RotationDefinition_BBOnly &def, double desiredAngle ); // set the torsion angle to the 'desiredAngle'

		std::vector<RotationDefinition_BBOnly> m_RotatePhi_BB; // Array: cache all the phi rotation atoms
		std::vector<RotationDefinition_BBOnly> m_RotatePsi_BB; // Array: cache all the psi rotation atoms
		std::vector<RotationDefinition_BBOnly> m_RotateOmega_BB; // Array: cache all the omega rotation atoms

	private:
		std::vector<Maths::dvector*> m_CAPos;
		std::vector<Maths::dvector*> m_BBPos;
		SegBreakType m_RotbreakType; ///< Are we rotating from the forward or backward anchor point?
	};


	//-------------------------------------------------
	/// \brief  BRIEF DESCRIPTION
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// Intended as a base class for both the ConformerBuilder
	/// and the SegPerturbation_BB. Used
	/// to apply backbone torsions to the particle system within
	/// the seg structure only, manipulating all atoms ...
	/// We cant use the rotatable bond list of the particle system,
	/// because it assumes that there are no seg breaks...
	/// \author  Jon Rea 
	class PD_API WholeBackboneTC
	{
	protected:

		//-------------------------------------------------
		/// \brief  BRIEF DESCRIPTION
		/// \details DETAILED USER'S DESCRIPTION
		///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
		/// \author  Jon Rea 
		class RotationDefinition_AllAtom: public RotationDefinition
		{
		public:
			RotationDefinition_AllAtom() :
			  RotationDefinition(),
				  conformerIndex(-1)
			  {
			  }

			  // extra variable required.
			  int conformerIndex;
		};

	public:
		WholeBackboneTC();

		void setBackBonePhiPsi( double _Phi, double _Psi );
		void setBackBoneOmega( double _Omg );

		void setMolecule( MoleculeBase& _MolBase, SegmentDef& _Range ); ///< Uses the ResidueRange derived SegmentDef which includes the boolean reverse term in isBackwards(). 
		void setMolecule( MoleculeBase& _MolBase, PickResidueRange& _Range, SegBreakType _breakType ); ///< Uses a standard ResidueRange and muct be told which anchor the segment is using.

		void validateAll() const; ///< Go through all definitions making sure they are valid

	protected:		
		void setBackBonePhiPsi_NoValidate( double _Phi, double _Psi ); ///< Dangerous call, only usable from derived classes that have performed checking of their own.

		std::vector<RotationDefinition_AllAtom> m_RotatePhi_AA; ///< cache all the phi rotation atoms
		std::vector<RotationDefinition_AllAtom> m_RotatePsi_AA; ///< cache all the psi rotation atoms
		std::vector<RotationDefinition_AllAtom> m_RotateOmega_AA; ///< cache all the omega rotation atoms

	private:
		void initialiseRotations(SegBreakType _breakType);		

		PickResidueRange m_ResRange;
		MoleculeBase* m_InnerMolBase;
	};


	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
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
	class PD_API SegWholeBackboneTC : public WholeBackboneTC, public virtual SegmentDefUser
	{
	public:
		SegWholeBackboneTC( SegmentDef& _segdef ) : 
		  SegmentDefUser(_segdef)
		  {
			  setMolecule(_segdef.getWorkSpace(),_segdef);
		  }
	};

	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// This class has to derive from SegmentDefUser and NOT,
	/// TCBase because it uses WorkSpaces RotBond array, not
	/// a custom internal array.
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API SegSidechainTC : public virtual SegmentDefUser
	{
	protected:

		struct RotBondCache // internal data structure
		{
		public:		
			RotBondCache();

			void calcCache( const WorkSpace& _WSpace, int sidechainResNum );
			inline size_t size() { return m_RotBondIndexes.size(); }
			inline int at(size_t bondIndex){ return m_RotBondIndexes[bondIndex]; } // return the rotbond index of bond 'bondIndex' in the residue

		private:
			inline void clear() { m_RotBondIndexes.clear(); }
			std::vector<int> m_RotBondIndexes; // array of indexes to Rotatable Bonds of a particlesystem
		};

	public:
		SegSidechainTC( SegmentDef& _segdef )
			: SegmentDefUser( _segdef )
		{
			initialiseRotations();
		}

	protected:
		void initialiseRotations();
		std::vector<RotBondCache> m_RotBonds; // a pointer array to the particle systems 'RotatableBond's that we want to use
	};


	//-------------------------------------------------
	//
	/// \brief  A class for applying desired torsions to the seg on demand.
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
#ifndef SWIG // Python really doesnt like diamond inheritance...
	class PD_API SegTorsionManipulator:
		public virtual SegmentDefUser,
		public SegSidechainTC,
		public SegWholeBackboneTC
	{
	public:
		SegTorsionManipulator( SegmentDef& _segdef );
	};
#endif

} // namespace 'Manipulator'

#endif

