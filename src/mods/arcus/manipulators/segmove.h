#ifndef LOOP_PERTURBATION_H
#define LOOP_PERTURBATION_H

#include <vector>

#include "manipulators/basicmoves.h"
#include "manipulators/rotamer_applicatorbase.h"
#include "workspace/segdef.h"
#include "segtorsions.h"
#include "pickers/pickbase.h"

namespace Manipulator
{
	//-------------------------------------------------
	//
	/// \brief  BRIEF DESCRIPTION
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// Interface Executive Base class for MonteCarlo-Type changes
	/// to current wspace structure of seg regionsonly; as opposed to
	/// StructurePertubation classes which is designed to function on the entire
	/// system.
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API SegPerturbation: public virtual SegmentDefUser, public MoveBase
	{
	public:
		SegPerturbation( SegmentDef &_segdef );
		virtual ~SegPerturbation();
		virtual SegPerturbation* clone() const = 0;

		virtual int apply() = 0; // (still pure virtual function) we are not implementing at this stage, wait for the derived class to implement it ...

	protected:
		Maths::FastRandom *m_Rand;
	};


	//-------------------------------------------------
	//
	/// \brief  Perturb the backbone torsions in a seg structure
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// Principle based on the class "Structure::Perturbation::BackboneTorsionalPertubation"
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API SegPerturbation_BBTorsional: public SegPerturbation, public SegWholeBackboneTC
	{
	public:
		SegPerturbation_BBTorsional( SegmentDef &_segdef, double probability, double maxMagnitude ):
		  SegmentDefUser( _segdef ),
			  SegPerturbation( _segdef ),
			  SegWholeBackboneTC( _segdef )
		  {
			  Set(probability, maxMagnitude);
		  }

		  virtual SegPerturbation_BBTorsional* clone() const { return new SegPerturbation_BBTorsional(*this); }

		  virtual int apply();

		  void Set( double probability, double maxMagnitude );
	private:
		double m_Probability; // the per-residue probability that a perturbation will occur
		double m_MaxMagnitude; // the magnitude of the perturbation
	};


	//-------------------------------------------------
	//
	/// \brief  Perturb the sidechain torsions in a seg structure
	///
	/// \details 
	/// In this class a RotBondCache is used to store the indexes of the relevent wspace rotbonds
	/// Principle based on the class "Structure::Perturbation::SidechainTorsionalPertubation"
	/// probability: The probability that the perturbation will occur
	/// magnitude: The magnitude of selected perturbations
	/// probability120: The probability that a large 120 degree perturbation will occur
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class PD_API SegPerturbation_SCTorsional: public SegPerturbation, public SegSidechainTC
	{
	public:
		SegPerturbation_SCTorsional( SegmentDef &_segdef, double probability, double maxMagnitude, double probability120 ):
		  SegmentDefUser( _segdef ),
			  SegPerturbation( _segdef ),
			  SegSidechainTC( _segdef )
		  {
			  Set( probability, maxMagnitude, probability120 ); // set the internal probabilities.
		  }

		  virtual SegPerturbation_SCTorsional* clone() const { return new SegPerturbation_SCTorsional(*this); }

		  virtual int apply();
		  void Set( double probability, double maxMagnitude, double probability120 ); // can dynamically set the internal probabilities

	private:
		double m_Probability; // the per-residue probability that a perturbation will occur
		double m_Probability120; // the per-residue probability that a large perturbation will occur
		double m_MaxMagnitude; // the magnitude of the perturbation
	};


	//-------------------------------------------------
	//
	/// \brief  Transform the sidechain atoms into a random rotamer conformation
	///
	/// \details DETAILED USER'S DESCRIPTION
	///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
	/// Principle based on the class "Structure::Perturbation::SidechainRotamerLibPertubation"
	/// probability: The probability that a rotamer perturbation will occur
	///
	/// \author  Jon Rea 
	///
	/// \todo STATE OF DEVELOPMENT
	///
	/// \bug BUGS?
	///
	class SegPerturbation_SCRotamer : public virtual SegmentDefUser, public RotamerApplicatorBase
	{
	public:
		SegPerturbation_SCRotamer( 
			SegmentDef &_segdef,
			const Library::RotamerLibrary &_rotlib,
			double probability, 
			double forcedChangeCutoff,
			RotamerMode _mode = ApplyCartesian );

		virtual SegPerturbation_SCRotamer* clone() const { return new SegPerturbation_SCRotamer(*this); }

		virtual int apply();

		void Set( double probability, double forcedChangeCutoff );

	protected:
		Maths::FastRandom *m_Rand;

	private:
		// private class members
		const Library::RotamerLibrary *rotlib;
		double m_Probability; // the per-residue probability that a rotamer replacement will occur
		double m_forcedchangecutoff;

		// private helper functions
		int getRotSterics(int ir, double *rotsterics, bool tryallrotamers);
		int estimateSCEntropy();
	};
}

#ifdef SWIG
%template(ObjectContainer_SegPerturbation) ObjectContainer<Manipulator::SegPerturbation>;
#endif

namespace Manipulator
{
	//-------------------------------------------------
	//
	/// \brief  The main class to hold and invoke the required perturbations
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
	class PD_API SegTweaker : public ObjectContainer< SegPerturbation >
	{
	public:
		SegTweaker();
		virtual SegTweaker* clone() const { return new SegTweaker(*this); }
		void applyPerturbations();
	};

} // namespace Manipulator

#endif

