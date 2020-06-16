#ifndef __BUDE_H
#define __BUDE_H

#include "forcefields/forcefield.h" // Provides the base class

class PD_API Particle;

namespace Physics
{
//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
/// Provides a place to store all of the Custom Properties from the Bude forcefield parameter file
/// The call to get the custom property of an atom
/// (eg. wspace.atom[i].getCustomProperty_double("HARDNESS");)
/// is slow, so it's best to get all the custom properties once in the forcefield setup() function
/// and store them all in a std::vector<BudeCustomProperties> rather than doing that within the inner loops.
///
/// \author Jon Crisp 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	struct BudeCustomProperties
	{
		BudeCustomProperties();

		char ElectrostaticType;       ///< N == 0, F == 1, D == 2, E == 3.
		double Hardness;              ///< Hardness of balls in steric calculation.
		double HydrophobicPotential;  ///< K, also known as Desolvation Potential. Used in the desolvation calculation.
		double DistNpNp;              ///< Cutoff for non-polar - non-polar interactions. Used in the desolvation calculation.
		double DistNpP;               ///< Cutoff for non-polar - polar interactions. Used in the desolvation calculation.
		double RadiusScaling;         ///< A factor for scaling of the atomic radius.
	};




//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Jon Crisp 
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class BudeForcefield: public ForcefieldBase
	{
	public:
		
		/// Constructor
		BudeForcefield( WorkSpace &newwspace );

		virtual BudeForcefield* clone() const
		{
			return new BudeForcefield(*this);
		}

		// This is a function you must implement - it will be called
		// automatically (by the ForcefieldContainer this forcefield component
		// is added to) once before anyone tries to use this forcefield
		// It should perform any essential setup functionality such as:
		//   - Reading in forcefield parameters
		//   - Setting up internal data strucutres (precalculated data,
		//       paramters, libraries, etc..)
		//   - Request a minimum cutoff from the neighbour list 
		// Important: This function must be able to be run multiple times, i.e.
		//  you must clear any data structures and recreate them rather than
		//  adding to them.
		// When this function finished
		virtual int setup();


		/// calcualtes energies AND forces
		/// Not used by BudeForcefield but has to be here since it's a virtual function from the Base class
		virtual void calcForces(){} 

		/// calcualtes energies only
		virtual void calcEnergies();   

		/// calculates and displays energies verbosely
		virtual void calcEnergiesVerbose(Verbosity::Type level); 
	
	// Set and Get functions
	public:

		///\brief Set the molecule index of the receptor within the workspace.
		/// eg. if the receptor is loaded first then it will have index 0.
		void setReceptorIndex(const int& _ReceptorIndex);
		const int& getReceptorIndex() const;

		///\brief Set the molecule index of the ligand within the workspace.
		/// eg. if the ligand is loaded second then it will have index 1.
		void setLigandIndex(const int& _LigandIndex);
		const int& getLigandIndex() const;

		void setCutoff(const double& _cutoff); ///< Set global interaction cutoff distance [Angstroms]
		const double& getCutoff() const;

		void setCutoff_elec_formal(const double& _cutoff_elec_formal); ///< Set cutoff for electrostatic calculations in the case of formal charge-charge interactions [Angstroms]
		const double& getCutoff_elec_formal() const;

		void setCutoff_elec_partial(const double& _cutoff_elec_partial); ///< Set cutoff for electrostatic calculations where the atom is defined as having a partial charge (Dipole interactions ie. H-bonding) [Angstroms]
		const double& getCutoff_elec_partial() const;

	protected:
		// Put private data here such as internal stores of current step number, state,
		// Library of structures, etc...
		//
		// e.g.:

		int m_ReceptorIndex;           ///< molecule index of receptor within WorkSpace
		int m_LigandIndex;             ///< molecule index of ligand within WorkSpace

		double m_Cutoff;               ///< Global interaction cutoff distance [Angstroms]

		// Would it be better to set these in the forcefield parameter file and have them as another atom property?
		double m_Cutoff_elec_formal;   ///< Cutoff for electrostatic calculations in the case of formal charge-charge interactions [Angstroms]
		double m_Cutoff_elec_partial;  ///< Cutoff for electrostatic calculations where the atom is defined as having a partial charge (Dipole interactions ie. H-bonding) [Angstroms]

		/// "Dielectric constant" (fudge factor) for electrostatic calculation.
		/// This is included based on previous parameterisation so that the
		/// electrostatic calculation should give a ~2.5kcal H-bond.
		/// Initialised in the constructor to (PhysicsConst::econv)/22.5.
		double m_Dielectric;

		/// Local store of the Custom Properties from the Bude forcefield parameter file.
		/// This way they can be accessed quickly by member functions within loops
		/// instead of having to call wspace.atom[i].getCustomProperty.
		std::vector<BudeCustomProperties> m_BudeCustomProperties;

		static const char m_elecType_N = 0;  ///< Type N will now be referred to as 0
		static const char m_elecType_F = 1;  ///< Type F will now be referred to as 1
		static const char m_elecType_D = 2;  ///< Type D will now be referred to as 2
		static const char m_elecType_E = 3;  ///< Type E will now be referred to as 3

		double m_EpotLigand_vdw;       ///< local store of LIGAND INTERNAL steric energy
		double m_EpotLigand_elec;      ///< local store of LIGAND INTERNAL electrostatic energy
		double m_EpotLigand_desolv;    ///< local store of LIGAND INTERNAL desolvation energy
		double m_EpotLigand_total;     ///< local store of LIGAND INTERNAL Total energy (Sum of steric, elec and desolv)

		double m_EpotComplex_vdw;      ///< local store of COMPLEX (Ligand-Receptor) steric energy
		double m_EpotComplex_elec;     ///< local store of COMPLEX (Ligand-Receptor) electrostatic energy
		double m_EpotComplex_desolv;   ///< local store of COMPLEX (Ligand-Receptor) desolvation energy
		double m_EpotComplex_total;    ///< local store of COMPLEX (Ligand-Receptor) Total energy (Sum of steric, elec and desolv)

		double m_Epot_vdw;             ///< local store of TOTAL steric energy
		double m_Epot_elec;            ///< local store of TOTAL electrostatic energy
		double m_Epot_desolv;          ///< local store of TOTAL desolvation energy
		double m_Epot_total;           ///< local store of TOTAL Total energy (Sum of steric, elec and desolv)

	protected:

		void resetLocalEnergies(); ///< Use this function to set the local energy stores to zero.

		void validateParams(const WorkSpace& _wspace) const;   ///< check that everything has been set properly by the user

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       
		
		/// prints the headers for the above function
		virtual void infoLineHeader() const; 

	};
}

#endif


