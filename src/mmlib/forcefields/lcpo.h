#ifndef __LCPO_H
#define __LCPO_H

#include "sasabase.h" // provides base class

// Maximal number of neighbors in the SASA Neighbor_list
#define Nlistmax 80

namespace Physics
{





//-------------------------------------------------
//
/// \brief Description: Surface Area calculations Using Numerical Counting, the linear
/// approximate LCPO algorithm [1] 
/// Energies/Forces are calculated using Atomic Solvation Paraemters
///
///
/// \details 
///    
/// References:
/// [1] Jorg Weiser, Peter S. Shenkin, W. Clark Still
/// Approximate Atomic Surfaces from Linear Combinations of Pairwise Overlaps (LCPO)
/// J. Comp. Chem., Vol. 20, No. 2, 217-230 (1999)
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
	class PD_API FF_SASA_LCPO : public ForcefieldBase 
		// NOTE: Should derive from public 'FF_SurfaceArea' when re-fitted
	{
	public:
		FF_SASA_LCPO(WorkSpace &newwspace) 
			: ForcefieldBase(newwspace)
		{
			name = "LCPO SASA";
			settodefault();
		}

		virtual FF_SASA_LCPO* clone() const { return new FF_SASA_LCPO(*this); }

		std::string ASPsection_name;
		double GlobalASP;
		int UpdateSasa;

		virtual void settodefault()
		{
			ASPsection_name = "";
			GlobalASP = 0.0;
			UpdateSasa = 5;
		};

	protected:
		virtual void setup();
		virtual void calcEnergiesVerbose(ForcefieldBase::AtomicVerbosity level);
		virtual void calcEnergies();
		virtual void calcForces();

		virtual void info() const;           ///< prints a little block of parameter information
		virtual void infoLine() const;       ///< prints a line of current energies
		virtual void infoLineHeader() const; ///< prints the headers for the above function






//-------------------------------------------------
//
/// \brief  BRIEF DESCRIPTION
///
/// \details DETAILED USER'S DESCRIPTION
///    DESCRIBE PURPOSE, INTERACTION WITH OTHER CLASSES, EXAMPLE CODE
///
/// \author Mike Tyka  
///
/// \todo STATE OF DEVELOPMENT
///
/// \bug BUGS?
///
		class PD_API SASA_Atom{
		public:		
			SASA_Atom():sigma(0){};
			int attype;               ///< link to the atom type
			double P1, P2, P3, P4;    ///< LCPO Parameters
			double radius;            ///< Vdw Radius + 1.4
			double SASA;              ///< SASA of individual atom
			Maths::dvector SASAderiv; ///< derivative of SASA (SASAderiv * sigma = force)
			double sigma;             ///< atomic solvation parameter
			signed char neighbors;    ///< Nr of neighbors |???
			char use;                 ///< include in calculation
		};






		class PD_API member
		{
		public:
			int size;
			int i[Nlistmax];         ///< index of neighbor
			double s[Nlistmax];      ///< this atoms overlap with its neighbor i
		};

	private:
		std::vector<SASA_Atom> SASAtype;
		std::vector<SASA_Atom> SASAatom;
		std::vector<member> Nlist;

		double epot_cav;
		double totalSASA;

		// can be used to implement forcefields like Eisenberg et al., Ooi et al etc..
		int readASPSection(const std::string &sectionname); // reads ASP parameters

		// this reads in the algorithm parameters (LCPO), called by setup();
		int readLCPOSection(); // reads the algorithm parameters

		// Solvent Accesable Surface Area calculations
		int calcNumericalSASA();
		int benchmarkNumericalSASA();

		int calcLCPOSASA();
		int calcLCPOSasaEnergies(bool recalc_nlist=true);
		int calcLCPOSasaForces_num();
		int calcLCPOSasaForces(bool dofullcalc);

		void testDerivatives();
	};


} // namespace Physics

#endif
