#ifndef __BUDE_EMC_H
#define __BUDE_EMC_H

#include <algorithm>
#include <vector>

#include "protocols/protocolbase.h" // Provides the base class

#include "typedefs.h"

#include "library/rotamerlib.h"
#include "manipulators/rotamer_applicatorbase.h"
#include "maths/fastrandom.h"
#include "workspace/workspace.fwd.h"

#include "bude_ff.h"

class PickResiduesFromFile;

namespace Protocol{

	/// Arrays of translations and rotations that each member of the Positional Descriptor represents
	struct PosDescTransRot
	{
		std::vector<double> transx;  ///< actual x translations
		std::vector<double> transy;  ///< actual y translations
		std::vector<double> transz;  ///< actual z translations
		std::vector<double> rotx;    ///< actual x rotations
		std::vector<double> roty;    ///< actual y rotations
		std::vector<double> rotz;    ///< actual z rotations
	};

	/// Positional Descriptor
	struct PosDesc
	{
		PosDesc();

		conformer_type tx;  ///< x translation
		conformer_type ty;  ///< y translation
		conformer_type tz;  ///< z translation
		conformer_type rx;  ///< x rotation
		conformer_type ry;  ///< y rotation
		conformer_type rz;  ///< z rotation

		//double score;

		//bool operator< (const PosDesc &rhs) const {return score < rhs.score;}
	};

	/// Store of information about the rotamer for a particular residue
	struct RotChoice
	{
		RotChoice();

		int ir; ///< sequential residue index (should this be residue number??)
		int irot; ///< rotamer number - refers to the relevant rotamer number in the chosen rotamer library (in BudeEMC).

		bool operator< (const RotChoice &rhs) const { return ir < rhs.ir; }
	};

	/// Conformer Descriptor - a container for a PosDesc, a ligand RotDesc and a receptor RotDesc.
	struct ConformerDescriptor
	{
		ConformerDescriptor();

		PosDesc PD; ///< Positional Descriptor
		std::vector<RotChoice> ligandRD; ///< Ligand Rotamer Descriptor
		std::vector<RotChoice> receptorRD; ///< Receptor Rotamer Descriptor

		double score;

		bool operator< (const ConformerDescriptor &rhs) const {return score < rhs.score;}
	};






//-------------------------------------------------
//
/// \brief  BUDE EMC protocol
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
	class PD_API BudeEMC : public ProtocolBase
	{
	public: 

		/// Protocols always take at least a pointer to a forcefield
		/// at creation time - that paramter is passed on to the 
		/// base class Protocol base, you dont need to worry about it
		BudeEMC( Physics::Forcefield& _ff, const Library::RotamerLibrary& _Lib, Manipulator::RotamerMode _mode = Manipulator::ApplyCartesian );
		virtual BudeEMC* clone() const 
		{ 
			return new BudeEMC(*this); 
		}

		/// Destructor - clean up your mess in here
		virtual ~BudeEMC();
		

		/// prints a little block of parameter information, the class
		/// should print its entire "state" through this call
		virtual void info() const; 
		
		
		/// this function is the main function which carries out the
		/// main loop of the simulation. It's outline implementation 
		/// is in example.cpp
		virtual int runcore(); 

	public:

		/// public data i.e. parameters. 
		
		/// Note that your Base class ProtcolBase already
		/// provides a number of public parameters which this class will inherit
		/// and that you can use directly (see protocol.h for more details)
		/// 
		/// int UpdateScr;   // screen/stdout update frequency (0=off)
		/// int UpdateTra;   // trajectory update frequency (0=off)
		/// int UpdateMon;   // monitor update frequency ( (0=off)
		/// int UpdateNList; // number of Steps between full neighborlist updates;
		/// Verbosity::Type OutputLevel;     // if true then completely silent mode (no screen output at all)
		/// size_t Steps;    // number of steps to be done

		int RandomSeed;           ///< random number seed

		/// Function for user to add possible translations and rotations to m_TransRotGrid.
		void setTransformations(const std::string& filename);

		/// Function for user to set the initial generation of PDs to use rather than using genInitialPDs().
		void setInitialGeneration(const std::string& filename);

		/// Set the location of the Rotamer Choices file for the ligand.
		/// Each residue listed in the file should be on a new line
		/// and each entry should contain (in this order):
		/// a chain ID, a residue number, a residue 3 letter code and a residue Icode (optional).
		void setLigandRotamerPerturbations(const std::string& filename);

		/// Set the location of the Rotamer Choices file for the receptor.
		/// Each residue listed in the file should be on a new line
		/// and each entry should contain (in this order):
		/// a chain ID, a residue number, a residue 3 letter code and a residue Icode (optional).
		void setReceptorRotamerPerturbations(const std::string& filename);

		/// Set the number of children that will initially be generated in the rotamer emc perturbation.
		void setInitialRotamerChildrenCount(size_t _initialRotamerChildrenCount);
		size_t getInitialRotamerChildrenCount() const; ///< Get the number of children that will initially be generated in the rotamer emc perturbation.

		/// Set the number of parents and children to use in the rotamer emc perturbation.
		void setRotamerGenerationSize(size_t _rotamerParentCount, size_t _rotamerChildrenCount);
		size_t getRotamerParentCount() const; ///< Get the number of parents used in the rotamer emc perturbation.
		size_t getRotamerChildrenCount() const; ///< Get the number of children used in the rotamer emc perturbation.
		size_t getRotamerGenerationSize() const; ///< Get the generation size used in the rotamer emc perturbation.

		/// Steps == Number of Generations
		/// use "Steps" rather than having a "size_t NumberOfGenerations"

		void setRotamerSteps(size_t _rotamerSteps); ///< Set the number of generations to use in the rotamer emc perturbation.
		size_t getRotamerSteps() const; ///< Get the number of generations that will be used in the rotamer emc perturbation.

		size_t InitialGenCount;   ///< number of initial parent Positional Descriptors to be generated
		size_t ChildrenPerParent; ///< number of children to generate per parent

		///\brief MutationRate determines how many of the Positional Descriptor's x, y and z translations
		/// and rotations are changed when generating a child from a parent.
		/// 
		/// eg. a mutation rate of 50% means that 3 out of 6 are varied ie. 777777 might become 547727.
		double MutationRate;  // NB should this be a size_t or a double?
		                      // Think it depends on whether want to say 0.5 or 50%.

		size_t NextGenCount;      ///< number of parents in all generations after the 1st

		int ReceptorIndex;     ///< molecule index of receptor within WorkSpace
		int LigandIndex;       ///< molecule index of ligand within WorkSpace

	protected:

		PosDescTransRot m_TransRotGrid; ///< store of translations and rotations represented by the values in a PosDesc.

		/// \brief A container for Conformer Descriptors
		/// This will allow Conformer Descriptors to be easily stored and then sorted
		std::vector<ConformerDescriptor> m_ConfDescs;

		/// Centre of geometry of ligand at initial pose
		Maths::dvector m_LigandCOG;

		/// Random number generator
		Maths::FastRandom m_Rand;

		/// Current step (generation) number
		size_t Step;

		/// Number of generations to use in the rotamer emc perturbation.
		size_t m_RotamerSteps;

		/// Current generation in the rotamer emc perturbation.
		size_t m_RotamerStep;

		/// Flag to keep track of whether the user has set the number of generations to do in the rotamer emc perturbation.
		bool m_DefaultRotamerSteps;

		/// Max values of integers used within the Positional Descriptors to represent translations and rotations.
		/// ie. number of possible x translations, y translations etc.
		size_t TxMax; ///< x Translations
		size_t TyMax; ///< y Translations
		size_t TzMax; ///< z Translations
		size_t RxMax; ///< x Rotations
		size_t RyMax; ///< y Rotations
		size_t RzMax; ///< z Rotations

		bool m_DoInitialGeneration; ///< False if the user has specified a list of initial PDs to use.

		std::vector<size_t> m_TotalPossRotamers; ///< Holds the number of possible rotamers for each residue in the WorkSpace.

		bool m_DefaultInitialRotamerChildrenCount; ///< Flag to keep track of whether the user has called setInitialRotamerChildrenCount().
		size_t m_InitialRotamerChildrenCount; ///< Number of children to generate in the initial step of the rotamer emc perturbation. Default == (NextGenCount*ChildrenPerParent).

		bool m_DefaultRotamerGenerationSize; ///< Flag to keep track of whether the user has called setRotamerGenerationSize().
		size_t m_RotamerParentCount;         ///< Number of parents to use in the rotamer emc perturbation.
		size_t m_RotamerChildrenCount;       ///< Number of children to use in the rotamer emc perturbation.
		size_t m_RotamerGenerationSize;      ///< Generation size in the rotamer emc perturbation (== m_RotamerParentCount + (m_RotamerParentCount * m_RotamerChildrenCount)).

		bool m_DoRotamerPerturbations; ///< Do EMC just with PosDescs or with rotamer perturbations too?
		bool m_DoLigandRotamers;
		bool m_DoReceptorRotamers;

		Manipulator::RandomRotamerApplicator m_RotlibApply; ///< Applies the rotamers
		std::string m_LigandRotChoiceFile;    ///< ligand choices file for rotamer perturbations.
		std::string m_ReceptorRotChoiceFile;    ///< receptor choices file for rotamer perturbations.

		/// Called by runcore() prior to execution validating the public parameters ReceptorIndex and LigandIndex.
		/// This allows checking that a receptor and a ligand actually exist within the WorkSpace.
		void validateParams(const WorkSpace& _wspace) const;

		///\brief Generate an initial list of Positional Descriptors using InitialGenCount then
		/// calculate their energies, sort them by that energy and only keep the top NextGenCount.
		void genInitialPDs(PosStore& _initLigandPose, WorkSpace& _wspace);

		///< Setup the residue pickers and all of the Rotamer Descriptors.
		void setupRotamers(WorkSpace& _wspace, PickResiduesFromFile& _ligandRotChoices, PickResiduesFromFile& _receptorRotChoices);

		void generateChildren();        ///< Generate children using ChildrenPerParent

		PosDesc mutatePD(const PosDesc& _PD); ///< Function to mutate parent PDs to generate children

		void mutatePDValue(const size_t _whichPDValue, PosDesc& _child); ///< Function to mutate the actual PD values eg. tx, ry etc.

		///\brief Apply the selected Positional Descriptor to the ligand's initial pose to generate the new pose
		/// for energy evaluation.
		///
		/// Have to apply the transformation (specified by the Positional Descriptor) to the initial pose when it
		/// is at the centre of geometry (COG) of the system. Therefore have to move m_LigandCOG to the COG of the
		/// system, then apply the transformation, then move it back to the result of the initial pose with the
		/// transformation applied.
		void make(PosStore& _posStore, const PosDesc& _PD);

		///\brief Rotamer EMC Perturbation function
		/// 
		/// Randomly mutate the ConformerDescriptor passed as an argument to generate a selection of
		/// child conformers. Score all of the children and choose the best conformer (ie. best ligandRD
		/// and receptorRD combination for this PD). Return this best conformer along with EITHER its
		/// score OR the average of the scores of the best N (eg. 10) conformers.
		ConformerDescriptor rotamerEmcPerturbation(
			WorkSpace& _wspace,
			PosStore& _initLigandPose,
			PickBase& _ligandPicker,
			PickBase& _receptorPicker,
			Protocol::ConformerDescriptor& _originalConfDesc );

		void generateInitialRDChildren(ConformerDescriptor& _originalConfDesc, std::vector<ConformerDescriptor>& _localConfDescChildren);

		void generateRDChildren(std::vector<ConformerDescriptor>& _localConfDescChildren);

		std::vector<RotChoice> mutateRotDesc(const std::vector<RotChoice>& _RD);

		void mutateRotChoice(const size_t _whichResidue, std::vector<RotChoice>& _child);

		/// Use overloaded "operator <" in ConformerDescriptor to sort the conformers by score.
		/// #include the STL <algorithm> file to allow use of the "sort" function.
		void sortConfDescs(std::vector<ConformerDescriptor>& _ConfDescs){sort(_ConfDescs.begin(),_ConfDescs.end());}

		/// Use overloaded "operator <" in RotChoice to sort the RDs by residue number (ir).
		void sortRD(std::vector<RotChoice>& _RD){sort(_RD.begin(),_RD.end());}

		///\brief Only keep NextGenCount poses/conformers for the next generation.
		void cullChildren(std::vector<ConformerDescriptor>& _ConfDescs,size_t _ParentCount){_ConfDescs.erase(_ConfDescs.begin() + _ParentCount,_ConfDescs.end());}

		/// Put private member function such as helpers and bits of algorithms
		/// etc.. here. You can use virtual functions to allow other people
		/// to modify your algorithm by overloading specific functional aspects 
		/// of your algorithm. Note that calling virtual functions does carry 
		/// a small overhead.

		/// prints a line of current energies/information/stepnumber etc..
		virtual void infoLine() const;       
		
		/// prints the headers for the above function
		virtual void infoLineHeader() const; 
	};
}

#endif



