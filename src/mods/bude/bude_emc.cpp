#include "global.h"

#include "manipulators/basicmoves.h"

#include "forcefields/forcefield.h"

#include "pickers/pickfromfile.h"

#include "workspace/neighbourlist.h"
#include "workspace/pospointer.h"
#include "workspace/workspace.h"

#include "bude_emc.h"

// Declare any namespaces you might need to use (optional)
//using namespace Physics;
using namespace Maths;
//using namespace Monitors;
// ...

// Protocols are all in one common namespace "Protocol"
namespace Protocol
{

	PosDesc::PosDesc()
	{
		tx = UCHAR_MAX;
		ty = UCHAR_MAX;
		tz = UCHAR_MAX;
		rx = UCHAR_MAX;
		ry = UCHAR_MAX;
		rz = UCHAR_MAX;
	}

	RotChoice::RotChoice()
	{
		ir = -1;
		irot = -1;
	}

	ConformerDescriptor::ConformerDescriptor()
	{
		score = -1.0;
	}

	BudeEMC::BudeEMC( Physics::Forcefield& _ff, const Library::RotamerLibrary& _Lib, Manipulator::RotamerMode _mode )
		: ProtocolBase(_ff),
		m_RotlibApply(_ff.getWSpace(),_Lib,_mode)
	{
		// You should also initilise your private/protected member variables
		RandomSeed = 894891;
		InitialGenCount = 0;
		ChildrenPerParent = 0;
		MutationRate = 0.8;
		NextGenCount = 0;

		m_DefaultInitialRotamerChildrenCount = true;
		m_InitialRotamerChildrenCount = 0;

		m_DefaultRotamerGenerationSize = true;
		m_RotamerParentCount = 0;
		m_RotamerChildrenCount = 0;
		m_RotamerGenerationSize = 0;

		m_DefaultRotamerSteps = true;
		m_RotamerSteps = 0;
		m_RotamerStep = 0;

		ReceptorIndex = -1;
		LigandIndex = -1;

		Step = 0;

		m_DoRotamerPerturbations = false;
		m_DoLigandRotamers = false;
		m_DoReceptorRotamers = false;

		// Set the default translations and rotations ( a 7x7x7x7x7x7 grid )

		m_TransRotGrid.transx.push_back(-3);
		m_TransRotGrid.transx.push_back(-2);
		m_TransRotGrid.transx.push_back(-1);
		m_TransRotGrid.transx.push_back(0);
		m_TransRotGrid.transx.push_back(1);
		m_TransRotGrid.transx.push_back(2);
		m_TransRotGrid.transx.push_back(3);
		m_TransRotGrid.rotx.push_back(DegToRad(-30));
		m_TransRotGrid.rotx.push_back(DegToRad(-20));
		m_TransRotGrid.rotx.push_back(DegToRad(-10));
		m_TransRotGrid.rotx.push_back(DegToRad(0));
		m_TransRotGrid.rotx.push_back(DegToRad(10));
		m_TransRotGrid.rotx.push_back(DegToRad(20));
		m_TransRotGrid.rotx.push_back(DegToRad(30));

		m_TransRotGrid.transy.push_back(-3);
		m_TransRotGrid.transy.push_back(-2);
		m_TransRotGrid.transy.push_back(-1);
		m_TransRotGrid.transy.push_back(0);
		m_TransRotGrid.transy.push_back(1);
		m_TransRotGrid.transy.push_back(2);
		m_TransRotGrid.transy.push_back(3);
		m_TransRotGrid.roty.push_back(DegToRad(-30));
		m_TransRotGrid.roty.push_back(DegToRad(-20));
		m_TransRotGrid.roty.push_back(DegToRad(-10));
		m_TransRotGrid.roty.push_back(DegToRad(0));
		m_TransRotGrid.roty.push_back(DegToRad(10));
		m_TransRotGrid.roty.push_back(DegToRad(20));
		m_TransRotGrid.roty.push_back(DegToRad(30));

		m_TransRotGrid.transz.push_back(-3);
		m_TransRotGrid.transz.push_back(-2);
		m_TransRotGrid.transz.push_back(-1);
		m_TransRotGrid.transz.push_back(0);
		m_TransRotGrid.transz.push_back(1);
		m_TransRotGrid.transz.push_back(2);
		m_TransRotGrid.transz.push_back(3);
		m_TransRotGrid.rotz.push_back(DegToRad(-30));
		m_TransRotGrid.rotz.push_back(DegToRad(-20));
		m_TransRotGrid.rotz.push_back(DegToRad(-10));
		m_TransRotGrid.rotz.push_back(DegToRad(0));
		m_TransRotGrid.rotz.push_back(DegToRad(10));
		m_TransRotGrid.rotz.push_back(DegToRad(20));
		m_TransRotGrid.rotz.push_back(DegToRad(30));

		//end

		TxMax = m_TransRotGrid.transx.size();
		TyMax = m_TransRotGrid.transy.size();
		TzMax = m_TransRotGrid.transz.size();
		RxMax = m_TransRotGrid.rotx.size();
		RyMax = m_TransRotGrid.roty.size();
		RzMax = m_TransRotGrid.rotz.size();

	}

	BudeEMC::~BudeEMC()
	{
		// If you need to create dynamically allocated memory, deallocate it here.
		// otherwise leave destructor blank.
	}


	void BudeEMC::info() const
	{
		printf("BudeEMC: \n");
	}

	void BudeEMC::setTransformations(const std::string& filename)
	{
		PosDescTransRot tempTransRotGrid;

		double valTransx = -1;
		double valTransy = -1;
		double valTransz = -1;
		double valRotx = -1;
		double valRoty = -1;
		double valRotz = -1;

		const char *whiteSpace = " \t\12\15"; ///< Whitespace can be "space", tab, newline character or end-of-line character
		const char *commentChar = "#\12\15";  ///< Possible comment characters

		std::ifstream transRotFile(filename.c_str(), std::ifstream::in);
		if( !transRotFile.is_open() ) throw IOException("Transformations file " + filename + " not found." );

		std::string line;

		while(std::getline(transRotFile, line))
		{
			removecomments(line, commentChar);
			std::vector<std::string> token;

			token = chopstr(line, whiteSpace);
			if(token.size() <= 0) continue; // this will be an empty line

			std::string transformation = token[0];

			if(cmpstring(transformation,"TRANSX"))
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valTransx);
					tempTransRotGrid.transx.push_back(valTransx);
				}
			else if(cmpstring(transformation,"TRANSY")) 
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valTransy);
					tempTransRotGrid.transy.push_back(valTransy);
				}
			else if(cmpstring(transformation,"TRANSZ")) 
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valTransz);
					tempTransRotGrid.transz.push_back(valTransz);
				}
			else if(cmpstring(transformation,"ROTX")) 
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valRotx);
					tempTransRotGrid.rotx.push_back(DegToRad(valRotx));
				}
			else if(cmpstring(transformation,"ROTY")) 
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valRoty);
					tempTransRotGrid.roty.push_back(DegToRad(valRoty));
				}
			else if(cmpstring(transformation,"ROTZ")) 
				for(size_t i = 1; i < token.size(); i++)
				{
					str2double(token[i],valRotz);
					tempTransRotGrid.rotz.push_back(DegToRad(valRotz));
				}
		}

		transRotFile.close();

		///\brief If tempTransRotGrid doesn't contain anything (ie. there is no user-defined transformations
		/// file) leave m_TransRotGrid set to default and throw an exception to warn the user.
		if( tempTransRotGrid.transx.size() > 0
			&& tempTransRotGrid.transy.size() > 0
			&& tempTransRotGrid.transz.size() > 0
			&& tempTransRotGrid.rotx.size() > 0
			&& tempTransRotGrid.roty.size() > 0
			&& tempTransRotGrid.rotz.size() > 0
			)
		{
			m_TransRotGrid = tempTransRotGrid;

			TxMax = m_TransRotGrid.transx.size();
			TyMax = m_TransRotGrid.transy.size();
			TzMax = m_TransRotGrid.transz.size();
			RxMax = m_TransRotGrid.rotx.size();
			RyMax = m_TransRotGrid.roty.size();
			RzMax = m_TransRotGrid.rotz.size();
		}


		else
			throw ParseException("Transformations File set incorrectly:\nthere must be at least one value set for each transformation.");
	}

	void BudeEMC::setInitialGeneration(const std::string& filename)
	{
		int xTrans = -1;
		int yTrans = -1;
		int zTrans = -1;
		int xRot = -1;
		int yRot = -1;
		int zRot = -1;

		const char *whiteSpace = " \t\12\15"; ///< Whitespace can be "space", tab, newline character or end-of-line character
		const char *commentChar = "#\12\15";  ///< Possible comment characters

		std::ifstream _initialPdChoices(filename.c_str(), std::ifstream::in);
		if( !_initialPdChoices.is_open() ) throw IOException("Transformations file " + filename + " not found." );

		std::string line;

		while(std::getline(_initialPdChoices, line))
		{
			removecomments(line, commentChar);
			std::vector<std::string> token;

			token = chopstr(line, whiteSpace);
			if(token.size() <= 0) continue; // this will be an empty line
			
			str2int(token[0],xTrans);
			str2int(token[1],yTrans);
			str2int(token[2],zTrans);
			str2int(token[3],xRot);
			str2int(token[4],yRot);
			str2int(token[5],zRot);

			ConformerDescriptor newconf;
			newconf.PD.tx = xTrans;
			newconf.PD.ty = yTrans;
			newconf.PD.tz = zTrans;
			newconf.PD.rx = xRot;
			newconf.PD.ry = yRot;
			newconf.PD.rz = zRot;
			m_ConfDescs.push_back(newconf);

		}

		_initialPdChoices.close();
		m_DoInitialGeneration = false;
	}

	void BudeEMC::setLigandRotamerPerturbations(const std::string& filename)
	{
		m_LigandRotChoiceFile = filename; // store the filename in a member variable.

		m_DoRotamerPerturbations = true; ///< if the user specifies a ligand rotamer choices file then they want to do some rotamer perturbations.
		m_DoLigandRotamers = true; ///< if the user specifies a ligand rotamer choices file then they want to do rotamer perturbations for the ligand.
	}

	void BudeEMC::setReceptorRotamerPerturbations(const std::string& filename)
	{
		m_ReceptorRotChoiceFile = filename; // store the filename in a member variable.

		m_DoRotamerPerturbations = true; ///< if the user specifies a receptor rotamer choices file then they want to do some rotamer perturbations.
		m_DoReceptorRotamers = true; ///< if the user specifies a receptor rotamer choices file then they want to do rotamer perturbations for the receptor.
	}

	void BudeEMC::setInitialRotamerChildrenCount(size_t _initialRotamerChildrenCount)
	{
		m_InitialRotamerChildrenCount = _initialRotamerChildrenCount;
		m_DefaultInitialRotamerChildrenCount = false;
	}

	size_t BudeEMC::getInitialRotamerChildrenCount() const
	{
		return m_InitialRotamerChildrenCount;
	}

	void BudeEMC::setRotamerGenerationSize(size_t _rotamerParentCount, size_t _rotamerChildrenCount)
	{
		m_RotamerParentCount = _rotamerParentCount;
		m_RotamerChildrenCount = _rotamerChildrenCount;

		m_RotamerGenerationSize = m_RotamerParentCount + (m_RotamerParentCount * m_RotamerChildrenCount);

		m_DefaultRotamerGenerationSize = false;

	}

	size_t BudeEMC::getRotamerChildrenCount() const
	{
		return m_RotamerParentCount;
	}

	size_t BudeEMC::getRotamerParentCount() const
	{
		return m_RotamerChildrenCount;
	}

	size_t BudeEMC::getRotamerGenerationSize() const
	{
		return m_RotamerGenerationSize;
	}

	void BudeEMC::setRotamerSteps(size_t _rotamerSteps)
	{
		m_RotamerSteps = _rotamerSteps;
		m_DefaultRotamerSteps = false;
	}

	size_t BudeEMC::getRotamerSteps() const
	{
		return m_RotamerSteps;
	}

	int BudeEMC::runcore()
	{
		printf("\nStarting BudeEMC protocol...\n");

		// Proxies

		WorkSpace& wspace = getWSpace();

		validateParams(wspace);

		// Set our random number seed
		m_Rand.reinitialise( RandomSeed );

		// Clear the store of number of possible rotamers for each residue in the WorkSpace
		m_TotalPossRotamers.clear();

		// Use a member posStore to store the initial ligand pose - this doesn't have to be a function call
		MoleculeRange& ligandMol = wspace.mol[LigandIndex];
		PickAtomRange ligandRange(wspace,ligandMol.ifirst,ligandMol.ilast);

		PosStore initLigandPose(wspace,ligandRange);
		initLigandPose.store(); ///< store the initial pose.

		m_LigandCOG = wspace.getCentreOfGeometry(ligandMol.ifirst,ligandMol.ilast);

		if( m_DoInitialGeneration )
		{
			genInitialPDs(initLigandPose,wspace); ///< push 'InitialGenCount' new PDs into m_PosDescs
		}

		/// Individual residue pickers for the ligand and receptor.
		/// Call the default constructor now so that the pickers are empty and
		/// then reset them below if the user wants to do rotamer perturbations.
		PickResiduesFromFile _ligandRotChoices(wspace);
		PickResiduesFromFile _receptorRotChoices(wspace);

		if( m_DoRotamerPerturbations )
		{
			setupRotamers(wspace, _ligandRotChoices, _receptorRotChoices);
		}

		/// Use "Step" to count generations, use Steps for the total number of generations
		for(Step = 0; Step < Steps; Step++)
		{
			printf("\n\nEMC Pose Search: Generation %d of %d \n",(Step+1),Steps); // Step+1 counts from eg. 1-20 rather than 0-19.

			if(Step > 0)
			{
				generateChildren(); // push yet more PDs into the m_PosDescs by randomly mutating each parent in turn
			}

			// Get the scores for each pose
			for(size_t i = 0; i < m_ConfDescs.size(); i++)
			{
				printf("\nPose %4d of %4d.\n",(i+1),m_ConfDescs.size()); // i+1 counts from eg. 1-30 rather than 0-29.

				// For the time being, don't calculate neighbour list at all
				// and assume that budeEMC will only ever be used with
				// bude forcefield.
				// refreshNeighborList(); --- dont do this, we ALWAYS want a full neighbour list generation because changes are large
				//wspace.nlist().calcNewList(); // Do this instead

				if( m_DoRotamerPerturbations ) // If we're doing rotamer perturbations
				{
					printf("\nDoing rotamer EMC perturbations:\n");

					ConformerDescriptor _tempConf;
					_tempConf = rotamerEmcPerturbation(wspace,initLigandPose,_ligandRotChoices,_receptorRotChoices,m_ConfDescs[i]);
					m_ConfDescs[i] = _tempConf;
				}
				else
				{
					if((i==0)&&(OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&(((Step) % UpdateScr) == 0))
					{
						//if(i==0)
						ff->infoLineHeader();
					}

					if(OutputLevel > Verbosity::Normal)
					{
						printf("Scoring pose %d without rotamer EMC perturbations...\n",(i+1));
					}

					make(initLigandPose, m_ConfDescs[i].PD); // Apply the current conformer
					ff->calcEnergies();
					m_ConfDescs[i].score = wspace.ene.epot;

					if((OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&(((Step) % UpdateScr) == 0))
					{
						ff->infoLine();
					}
				}
				
				runmonitors();
			}
			
			//// Display the infoline every so often (every UpdateScr generations)
			//if((OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&
			//	 (((Step) % UpdateScr) == 0))
			//{
			//	infoLine();
			//}

			printf("\nSorting poses and culling unwanted children to keep best %2d...\n",NextGenCount);

			sortConfDescs(m_ConfDescs); // Sort the m_ConfDescs using a pre-defined standard function call in #include <algorithm>

			// Check that there are currently more than NextGenCount Conformer Descriptors in m_ConfDescs before culling children.
			if( m_ConfDescs.size() > NextGenCount )
			{
				cullChildren(m_ConfDescs,NextGenCount); // Remove conformers from NextGenCount to the end of the array
			}
			else if( m_ConfDescs.size() < NextGenCount )
			{
				printf("\n\tWARNING:\n");
				printf("\t\tThe number of poses (%d) being carried on to the next generation (Gen %d of %d) is less than NextGenCount (%d).\n",m_ConfDescs.size(),(Step+1),Steps,NextGenCount);
				printf("\t\tIs this what the user desires?\n\n");
			}

			if(OutputLevel > Verbosity::Silent)
			{
				printf("Total Potential Energy of chosen poses (best first):\n");
				for(size_t i = 0; i < m_ConfDescs.size(); i++)
				{
					printf("\tPose %4d:  %10.4f",(i+1),m_ConfDescs[i].score);
					printf(" \tPD: %d%d%d%d%d%d",m_ConfDescs[i].PD.tx,m_ConfDescs[i].PD.ty,m_ConfDescs[i].PD.tz,m_ConfDescs[i].PD.rx,m_ConfDescs[i].PD.ry,m_ConfDescs[i].PD.rz);
					printf(" \tLRD:");
					for(size_t j = 0; j < m_ConfDescs[i].ligandRD.size(); j++)
					{
						printf(" %2d",m_ConfDescs[i].ligandRD[j].irot);
					}
					printf(" \tRRD:");
					for(size_t j = 0; j < m_ConfDescs[i].receptorRD.size(); j++)
					{
						printf(" %2d",m_ConfDescs[i].receptorRD[j].irot);
					}
					printf("\n");
				}
			}

			printf("\nAppending best pose to trajectory...\n");

			//// Add the chosen parent structures to the output trajectory
			//for(size_t i = 0; i < m_ConfDescs.size(); i++)

			// Apply only the best pose/conformer
			// NB This also means that at the end of the BudeEMC protocol the WorkSpace will
			// contain this best pose/conformer which may be useful so that the user can
			// immediately start a new BudeEMC protocol with a finer grid for example.
			for(size_t i = 0; i < 1; i++)
			{
				make(initLigandPose, m_ConfDescs[i].PD); // Apply the current conformer

				// Apply the current LIGAND rotamer descriptor
				for(size_t j = 0; j < m_ConfDescs[i].ligandRD.size(); j++)
				{
					int _rotamerNumber = m_ConfDescs[i].ligandRD[j].irot;

					// Make sure a rotamer has been set for this residue
					if( _rotamerNumber != -1 )
					{
						// old version -> m_Rotlib.changeSidechainConformation(&wspace,m_ConfDescs[i].ligandRD[j].ir,_rotamerNumber);
						m_RotlibApply.applyRotamer( m_ConfDescs[i].ligandRD[j].ir,_rotamerNumber );
					}
				}

				// Apply the current RECEPTOR rotamer descriptor
				for(size_t j = 0; j < m_ConfDescs[i].receptorRD.size(); j++)
				{
					int _rotamerNumber = m_ConfDescs[i].receptorRD[j].irot;
					
					// Make sure a rotamer has been set for this residue
					if( _rotamerNumber != -1 )
					{
						// old version -> m_Rotlib.changeSidechainConformation(&wspace,m_ConfDescs[i].receptorRD[j].ir,_rotamerNumber);
						m_RotlibApply.applyRotamer( m_ConfDescs[i].receptorRD[j].ir, _rotamerNumber );
					}
				}

				wspace.outtra.append();
			}
		}

		printf("\nEnd of BudeEMC protocol.\n");

		return Step;
	}

	void BudeEMC::validateParams(const WorkSpace& _wspace) const
	{
		if( (LigandIndex >= 0) && (ReceptorIndex >= 0)
			&& (LigandIndex != ReceptorIndex)
			&& (LigandIndex < _wspace.nMolecules()) && (ReceptorIndex < _wspace.nMolecules()) )
			return;

		else
		{
			throw ProcedureException("Ligand Index and Receptor Index have been set incorrectly");
		}
	}

	void BudeEMC::genInitialPDs(PosStore& _initLigandPose, WorkSpace& _wspace)
	{
		// In case runcore is called multiple times, clear m_ConfDescs so we start with a clean slate.
		m_ConfDescs.clear();

		printf("\nGenerating initial %d Positional Descriptors...\n",InitialGenCount);

		// If the user hasn't set InitialGenCount then warn them since otherwise no initial
		// parent PDs will be generated which will mean that no PDs will ever be generated.
		if( InitialGenCount == 0 )
		{
			throw ProcedureException("\nInitialGenCount has not been set, or has been set to 0. This will prevent any EMC pose search from taking place.\n");
		}

		for(size_t i = 0; i < InitialGenCount; i++)
		{
			PosDesc newpos;
			ConformerDescriptor newconf;
			newconf.PD.tx = m_Rand.next(TxMax);
			newconf.PD.ty = m_Rand.next(TyMax);
			newconf.PD.tz = m_Rand.next(TzMax);
			newconf.PD.rx = m_Rand.next(RxMax);
			newconf.PD.ry = m_Rand.next(RyMax);
			newconf.PD.rz = m_Rand.next(RzMax);
			m_ConfDescs.push_back(newconf);
		}
	}

	void BudeEMC::setupRotamers(WorkSpace& _wspace, PickResiduesFromFile& _ligandRotChoices, PickResiduesFromFile& _receptorRotChoices)
	{
		// First get the number of possible rotamers for each residue in the WorkSpace.
		// Do this outside the individual ligand/receptor loops because residues must
		// be in exactly the same order as they are in the workspace.
		for( size_t ir = 0; ir < _wspace.nResidues(); ir++ )
		{
			// old -> m_TotalPossRotamers.push_back( m_Rotlib.getRotamerCount(_wspace.getAtom(_wspace.getRes(ir).ifirst).parentl3name) );
			m_TotalPossRotamers.push_back( m_RotlibApply.nRot( ir ) );
		}

		// Set up local Rotamer Descriptors since at this stage every Conformer Descriptor will be given
		// the same ligandRD and receptorRD which will contain residue indices but not rotamer numbers.
		std::vector<RotChoice> _tempLigandRD;
		_tempLigandRD.clear();   // clear it in case runcore() is called multiple times
		std::vector<RotChoice> _tempReceptorRD;
		_tempReceptorRD.clear(); // clear it in case runcore() is called multiple times

		if( m_DoLigandRotamers )
		{
			/// If user wants to do rotamer perturbations for the ligand
			/// reset the picker with the appropriate filename.
			_ligandRotChoices = PickResiduesFromFile(_wspace, m_LigandRotChoiceFile);

			const int _ligandStartRes = _wspace.mol[LigandIndex].irfirst; ///< local store of the INDEX of the first ligand residue.
			const int _ligandEndRes = _wspace.mol[LigandIndex].irlast;    ///< local store of the INDEX of the last ligand residue.

			// loop over all ligand residues
			for(int ir = _ligandStartRes; ir <= _ligandEndRes; ir++)
			{
				// if this residue has been picked for variation by the user
				if( _ligandRotChoices.matches( _wspace.res[ir] ) )
				{
					RotChoice _tempRotChoice;
					_tempRotChoice.ir = ir; // make a RotChoice for this residue

					_tempLigandRD.push_back( _tempRotChoice ); // store the RotChoice in a Rotamer Descriptor
				}
			}

			// not sure if this is strictly necessary, but sort so the RotChoices in the Rotamer Descriptor are ordered by residue index.
			sortRD(_tempLigandRD);
		}

		if( m_DoReceptorRotamers )
		{
			/// If user wants to do rotamer perturbations for the receptor
			/// reset the picker with the appropriate filename.
			_receptorRotChoices = PickResiduesFromFile(_wspace, m_ReceptorRotChoiceFile);

			const int _receptorStartRes = _wspace.mol[ReceptorIndex].irfirst; ///< local store of the INDEX of the first receptor residue.
			const int _receptorEndRes = _wspace.mol[ReceptorIndex].irlast;    ///< local store of the INDEX of the last receptor residue.

			// loop over all receptor residues
			for(int ir = _receptorStartRes; ir <= _receptorEndRes; ir++)
			{
				// if this residue has been picked for variation by the user
				if( _receptorRotChoices.matches( _wspace.res[ir] ) )
				{
					RotChoice _tempRotChoice;
					_tempRotChoice.ir = ir; // make a RotChoice for this residue

					_tempReceptorRD.push_back( _tempRotChoice ); // store the RotChoice in a Rotamer Descriptor
				}
			}

			// not sure if this is strictly necessary, but sort so the RotChoices in the Rotamer Descriptor are ordered by residue index.
			sortRD(_tempReceptorRD);
		}

		// Fill up the ligandRD and receptorRD of every Conformer Descriptor generated by genInitialPDs()
		// according to the residues that the user has chosen for variation.
		// The number of "RotChoice" in each RD will be equal to the number of residues chosen for ligand
		// or receptor. Each RotChoice will have the correct residue index (ir) but the rotamer number (irot)
		// will not have been set yet (and will therefore be -1).
		for(size_t i = 0; i < m_ConfDescs.size(); i++)
		{
			m_ConfDescs[i].ligandRD = _tempLigandRD;
			m_ConfDescs[i].receptorRD = _tempReceptorRD;
		}

		//
		// Finally check whether the user has set EMC parameters and if not set them to the default:
		//

		// The default number of generations (Steps) is the same as for the PosDesc EMC procedure.
		if(m_DefaultRotamerSteps)
		{
			m_RotamerSteps = Steps;
		}

		// The default for m_InitialRotamerChildrenCount is the same number of children
		// generated for the PosDesc EMC procedure (NextGenCount*ChildrenPerParent).
		if(m_DefaultInitialRotamerChildrenCount)
		{
			m_InitialRotamerChildrenCount = NextGenCount*ChildrenPerParent;
		}

		// The default for m_RotamerParentCount is the number of parents used in the PosDesc EMC procedure (NextGenCount).
		// The default for m_RotamerChildrenCount is the number of children used in the PosDesc EMC procedure (ChildrenPerParent).
		// m_RotamerGenerationSize is therefore the generation size in the PosDesc EMC procedure (== NextGenCount + (NextGenCount * ChildrenPerParent)).
		if(m_DefaultRotamerGenerationSize)
		{
			m_RotamerParentCount = NextGenCount;
			m_RotamerChildrenCount = ChildrenPerParent;
			m_RotamerGenerationSize = m_RotamerParentCount + (m_RotamerParentCount * m_RotamerChildrenCount);
		}
	}

	void BudeEMC::generateChildren()
	{
		printf("\nGenerating PD children.....\n");

		// Can't just use NextGenCount here in case there aren't that many conformers in m_ConfDescs.
		// This may happen when Step == 1 if (the number of PDs specified by setInitialGeneration() + 1)
		// is less than NextGenCount.
		//
		// eg. if the number of PDs specified by setInitialGeneration() is 13, m_ConfDescs.size() will be 14 at
		// this point so if NextGenCount is 15 then the first child (m_ConfDescs[14]) will be treated as the
		// 15th parent and used to generate children itself. This would be bad!
		size_t _numberOfParents = m_ConfDescs.size();

		for(size_t i = 0; i < _numberOfParents; i++)
		{
			if(OutputLevel > Verbosity::Normal)
			{
				printf("Parent number: %d \n",(i+1));
			}

			for(size_t j = 0; j < ChildrenPerParent; j++)
			{
				if(OutputLevel > Verbosity::Normal)
				{
					printf("\tChild number: %d \n",(j+1));
				}

				ConformerDescriptor newconf;

				// Must keep the ligandRD and receptorRD of the parent conformer when mutating the PD.
				newconf.ligandRD = m_ConfDescs[i].ligandRD;
				newconf.receptorRD = m_ConfDescs[i].receptorRD;

				newconf.PD = mutatePD(m_ConfDescs[i].PD);
				m_ConfDescs.push_back( newconf );
			}
		}

		printf("....done\n");
	}

	PosDesc BudeEMC::mutatePD(const PosDesc& _PD)
	{
		/// Child Positional Descriptor produced by mutating the parent passed as an argument.
		/// Initialise as the parent PosDesc
		PosDesc child = _PD;

		size_t NumberOfTransformations = 6; ///< x, y and z translations and rotations (3 translations, 3 rotations)

		///\brief Number of mutations made will be from 1 up to the percentage specified by MutationRate
		size_t NumberOfMutations = (size_t)m_Rand.next(1,((size_t)(MutationRate*(double)NumberOfTransformations)+1));
		
		int mutated = 0; ///< Used as a counter to keep track of number of PD values mutated

		std::vector<bool> taken; ///< Keep track of which PD values have already been mutated
		taken.resize(NumberOfTransformations,false);   ///< Initialise all to false

		D_ASSERT(NumberOfMutations < taken.size(), 
				CodeException, 
				"There is an error in the code, the condition 'NumberOfMutations < taken.size()' must be true.");

		while( mutated < NumberOfMutations )
		{
			size_t WhichPDValue = (size_t)m_Rand.next(NumberOfTransformations);

			if( ! taken[WhichPDValue] )
			{
				mutatePDValue(WhichPDValue,child);
				taken[WhichPDValue] = true;
				mutated++;
			}
		}
		return child;
	}

	void BudeEMC::mutatePDValue(const size_t _whichPDValue, PosDesc& _child)
	{
		switch(_whichPDValue)
		{
			case 0:
				_child.tx = m_Rand.next(TxMax);
				break;
			case 1:
				_child.ty = m_Rand.next(TyMax);
				break;
			case 2:
				_child.tz = m_Rand.next(TzMax);
				break;
			case 3:
				_child.rx = m_Rand.next(RxMax);
				break;
			case 4:
				_child.ry = m_Rand.next(RyMax);
				break;
			case 5:
				_child.rz = m_Rand.next(RzMax);
				break;
			default:
				THROW(CodeException,"Chosen PD value for mutation is out of range. Only 6 values can currently exist");

		}
	}

	void BudeEMC::make(PosStore& _posStore, const PosDesc& _PD)
	{
		_posStore.revert();
		
		//printf("Translations:\nx:%d y:%d z:%d\n",_PD.tx,_PD.ty,_PD.tz);

		dvector move(m_TransRotGrid.transx[_PD.tx], m_TransRotGrid.transy[_PD.ty], m_TransRotGrid.transz[_PD.tz]);
		//move.setTo(0,0,0);
		matrix3x3 mat;

		//printf("Rotations:\nx:%d y:%d z:%d\n",_PD.rx,_PD.ry,_PD.rz);

		mat.setToXYZrot(m_TransRotGrid.rotx[_PD.rx], m_TransRotGrid.roty[_PD.ry], m_TransRotGrid.rotz[_PD.rz]);

		move.add( m_LigandCOG );
		for( size_t i = 0; i < _posStore.size(); i++ )
		{
			_posStore.p(i).sub( m_LigandCOG );
			_posStore.p(i).mulmat(mat);
			_posStore.p(i).add(move);			
		}
	}

	ConformerDescriptor BudeEMC::rotamerEmcPerturbation(
		WorkSpace& _wspace,
		PosStore& _initLigandPose,
		PickBase& _ligandPicker,
		PickBase& _receptorPicker,
		Protocol::ConformerDescriptor& _originalConfDesc )
	{
		std::vector<ConformerDescriptor> _localConfDescChildren;

		_localConfDescChildren.push_back(_originalConfDesc); // keep the original as well as mutating it to generate children.

		generateInitialRDChildren(_originalConfDesc,_localConfDescChildren);

		/// Use "Step" to count generations, use Steps for the total number of generations
		for(m_RotamerStep = 0; m_RotamerStep < m_RotamerSteps; m_RotamerStep++)
		{
			printf("\n\nEMC Conformer Search: Generation %d of %d \n",(m_RotamerStep+1),m_RotamerSteps); // m_RotamerStep+1 counts from eg. 1-20 rather than 0-19.

			if( m_RotamerStep > 0 )
			{
				generateRDChildren(_localConfDescChildren);
			}

			// Apply each conformer and do energy calculation.
			for(size_t i = 0; i < _localConfDescChildren.size(); i++)
			{
				if( i==0 )
				{
					printf("\nApplying Rotamer Descriptors and calculating energies for %d conformers...\n",_localConfDescChildren.size());
					if((OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&(((Step) % UpdateScr) == 0))
					{
						ff->infoLineHeader();
					}
				}
				
				if(OutputLevel > Verbosity::Normal)
				{
					printf("Conformer: %d \n",(i+1)); // i+1 counts from eg. 1-30 rather than 0-29.
				}

				// Apply the current positional descriptor
				make(_initLigandPose, _localConfDescChildren[i].PD);

				// Apply the current LIGAND rotamer descriptor
				for(size_t j = 0; j < _localConfDescChildren[i].ligandRD.size(); j++)
				{
					int _rotamerNumber = _localConfDescChildren[i].ligandRD[j].irot;

					// Make sure a rotamer has been set for this residue
					if( _rotamerNumber != -1 )
					{
						// old version -> m_Rotlib.changeSidechainConformation(&_wspace,_localConfDescChildren[i].ligandRD[j].ir,_rotamerNumber);
						m_RotlibApply.applyRotamer( _localConfDescChildren[i].ligandRD[j].ir, _rotamerNumber );
					}
				}

				// Apply the current RECEPTOR rotamer descriptor
				for(size_t j = 0; j < _localConfDescChildren[i].receptorRD.size(); j++)
				{
					int _rotamerNumber = _localConfDescChildren[i].receptorRD[j].irot;
					
					// Make sure a rotamer has been set for this residue
					if( _rotamerNumber != -1 )
					{
						// old version -> m_Rotlib.changeSidechainConformation(&_wspace,_localConfDescChildren[i].receptorRD[j].ir,_rotamerNumber);
						m_RotlibApply.applyRotamer( _localConfDescChildren[i].receptorRD[j].ir, _rotamerNumber );
					}
				}

				ff->calcEnergies();
				_localConfDescChildren[i].score = _wspace.ene.epot;

				if((OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&(((Step) % UpdateScr) == 0))
					{
						ff->infoLine();
					}
			}

			printf("\n....done\n\nSorting conformers and culling unwanted children to keep best %2d...\n",m_RotamerParentCount);

			sortConfDescs(_localConfDescChildren);

			// Make sure that _localConfDescChildren actually contains more than
			// m_RotamerParentCount conformers before trying to cull the rest.
			// If there are exactly m_RotamerParentCount then there's no need to remove any.
			if( _localConfDescChildren.size() > m_RotamerParentCount )
			{
				cullChildren(_localConfDescChildren,m_RotamerParentCount);
			}
			else if( _localConfDescChildren.size() < m_RotamerParentCount )
			{
				printf("\n\tWARNING:\n");
				printf("\t\tThe number of conformers (%d) being carried on to the next generation (Gen %d of %d) is less than the RotamerParentCount (%d).\n",_localConfDescChildren.size(),(m_RotamerStep+1),m_RotamerSteps,m_RotamerParentCount);
				printf("\t\tIs this what the user desires?\n\n");
			}

			printf("....done\n\n");

			if(OutputLevel > Verbosity::Quiet)
			{
				printf("Total Potential Energy of chosen conformers (best first):\n");
				for(size_t i = 0; i < _localConfDescChildren.size(); i++)
				{
					printf("\tConformer %4d:  %10.4f",(i+1),_localConfDescChildren[i].score);
					printf(" \tPD: %d%d%d%d%d%d",_localConfDescChildren[i].PD.tx,_localConfDescChildren[i].PD.ty,_localConfDescChildren[i].PD.tz,_localConfDescChildren[i].PD.rx,_localConfDescChildren[i].PD.ry,_localConfDescChildren[i].PD.rz);
					printf(" \tLRD:");
					for(size_t j = 0; j < _localConfDescChildren[i].ligandRD.size(); j++)
					{
						printf(" %2d",_localConfDescChildren[i].ligandRD[j].irot);
					}
					printf(" \tRRD:");
					for(size_t j = 0; j < _localConfDescChildren[i].receptorRD.size(); j++)
					{
						printf(" %2d",_localConfDescChildren[i].receptorRD[j].irot);
					}
					printf("\n");
				}
			}
		}

		if(OutputLevel > Verbosity::Silent)
		{
			printf("\nPicking best conformer after rotamer emc perturbations.\n");

			printf(" \tPD: %d%d%d%d%d%d",_localConfDescChildren[0].PD.tx,_localConfDescChildren[0].PD.ty,_localConfDescChildren[0].PD.tz,_localConfDescChildren[0].PD.rx,_localConfDescChildren[0].PD.ry,_localConfDescChildren[0].PD.rz);
			printf(" \tLRD:");
			for(size_t j = 0; j < _localConfDescChildren[0].ligandRD.size(); j++)
			{
				printf(" %2d",_localConfDescChildren[0].ligandRD[j].irot);
			}
			printf(" \tRRD:");
			for(size_t j = 0; j < _localConfDescChildren[0].receptorRD.size(); j++)
			{
				printf(" %2d",_localConfDescChildren[0].receptorRD[j].irot);
			}
			printf("\n");
		}

		return _localConfDescChildren[0]; // returns the ConformerDescriptor with the lowest energy
	}

	void BudeEMC::generateInitialRDChildren(ConformerDescriptor& _originalConfDesc, std::vector<ConformerDescriptor>& _localConfDescChildren)
	{
		printf("\nGenerating initial rotamer children....\n");

		for(size_t i = 0; i < m_InitialRotamerChildrenCount; i++)
		{
			if(OutputLevel > Verbosity::Normal)
			{
				printf("\tChild number: %d \n",(i+1));
			}

			ConformerDescriptor newconf; // a local ConformerDescriptor

			if( m_DoLigandRotamers )
			{
				newconf.ligandRD = mutateRotDesc(_originalConfDesc.ligandRD); // mutate the ligand rotamers of this conformer
			}
			if( m_DoReceptorRotamers )
			{
				newconf.receptorRD = mutateRotDesc(_originalConfDesc.receptorRD); // mutate the receptor rotamers of this conformer
			}

			_localConfDescChildren.push_back( newconf ); // add the new child to the store of rotamer children
		}

		// Give each of the children the PD from the original ConformerDescriptor
		for(size_t i = 0; i < _localConfDescChildren.size(); i++)
		{
			_localConfDescChildren[i].PD = _originalConfDesc.PD;
		}

		printf("....done\n");
	}

	void BudeEMC::generateRDChildren(std::vector<ConformerDescriptor>& _localConfDescChildren)
	{
		printf("\nGenerating rotamer children....\n");

		// Can't just use m_RotamerParentCount here in case there aren't that many conformers in _localConfDescChildren.
		// This may happen when m_RotamerStep == 1 if (m_InitialRotamerChildrenCount + 1) < m_RotamerParentCount.
		//
		// eg. if m_InitialRotamerChildrenCount is 10, _localConfDescChildren.size() will be 11 at this point
		// so if m_RotamerParentCount is 12 then the first child (_localConfDescChildren[11]) will be treated as the 12th
		// parent and used to generate children itself. This would be bad!
		size_t _numberOfParents = _localConfDescChildren.size();

		// Generate m_RotamerChildrenCount children for each parent already in _localConfDescChildren.
		for(size_t i = 0; i < _numberOfParents; i++)
		{
			if(OutputLevel > Verbosity::Normal)
			{
				printf("\tParent number: %d \n",(i+1));
			}

			for(size_t j = 0; j < m_RotamerChildrenCount; j++)
			{
				if(OutputLevel > Verbosity::Normal)
				{
					printf("\tChild number: %d \n",(j+1));
				}

				ConformerDescriptor newconf; // a local ConformerDescriptor

				// Must keep the PD of the parent
				newconf.PD = _localConfDescChildren[i].PD;

				if( m_DoLigandRotamers )
				{
					newconf.ligandRD = mutateRotDesc(_localConfDescChildren[i].ligandRD); // mutate the ligand rotamers of this conformer
				}
				if( m_DoReceptorRotamers )
				{
					newconf.receptorRD = mutateRotDesc(_localConfDescChildren[i].receptorRD); // mutate the receptor rotamers of this conformer
				}

				_localConfDescChildren.push_back( newconf ); // add the new child to the store of rotamer children
			}
		}
		printf("....done\n");
	}

	std::vector<RotChoice> BudeEMC::mutateRotDesc(const std::vector<RotChoice>& _RD)
	{
		/// Child Rotamer Descriptor that will be produced by mutating the parent passed as an argument.
		/// Initialise as the parent Rotamer Descriptor.
		std::vector<RotChoice> child = _RD;

		size_t NumberOfResidues = child.size(); ///< Number of residues that the user has said can be mutated.

		/// Number of mutations made will be from 1 up to the percentage specified by MutationRate
		size_t NumberOfMutations = (size_t)m_Rand.next(1,((size_t)(MutationRate*(double)NumberOfResidues)+1));
		
		int mutated = 0; ///< Used as a counter to keep track of number of PD values mutated

		std::vector<bool> taken; ///< Keep track of which PD values have already been mutated
		taken.resize(NumberOfResidues,false);   ///< Initialise all to false

		D_ASSERT(NumberOfMutations < taken.size(), 
				CodeException, 
				"There is an error in the code, the condition 'NumberOfMutations < taken.size()' must be true.");

		while( mutated < NumberOfMutations )
		{
			size_t WhichResidue = (size_t)m_Rand.next(NumberOfResidues);

			if( ! taken[WhichResidue] )
			{
				mutateRotChoice(WhichResidue,child);
				taken[WhichResidue] = true;
				mutated++;
			}
		}
		return child;
	}

	void BudeEMC::mutateRotChoice(const size_t _whichResidue, std::vector<RotChoice>& _child)
	{
		// Randomly mutate the rotamer number of the chosen residue to a value
		// between 0 and the total number of rotamers in the library for that residue.
		_child[_whichResidue].irot = m_Rand.next( m_TotalPossRotamers[_child[_whichResidue].ir] );
	}

		//// Main simulation loop, e.g. loops of the number of steps
		//// to be simulated. The protocol here is a dummy protocol and serves
		//// exemplary function only.
		//for(Step = 0; Step < Steps; Step++) {

		//	// refresh the neighbor list - this function is from the base class
		//	// ProtocolBase and will make sure the neighbor list in wspace
		//	// is uptodate according to the update frequence (UpdateNList)
		//	// that the user has set 
		//	refreshNeighborList();

		//	// Calculate the forces or energies of the system
		//	ff->calcForces();
		//	// If you dont need the forces but just the energies then just
		//	// call ff->calcEnergies(); instead, that is sometimes faster

		//	// The total potential energy will we available in
		//	// getWSpace().ene.epot

		//	// Call any monitors that might have been added to you as a protocol
		//	runmonitors();

		//	// Now apply the forces/change the system/etc..
		//	// The current atom positions  are in getWSpace().cur.atom[AtomNumber].p.{x/y/z}
		//	// The current atom forces     are in getWSpace().cur.atom[AtomNumber].f.{x/y/z}
		//	// The current atom velocities are in getWSpace().cur.atom[AtomNumber].v.{x/y/z}
		//	// getWSpace().old.atom[AtomNumber].{p/f/v} can be used to store these quantities
		//	// from the last time step if required. Alternatively a local isntance of
		//	// SnapShot can be used for the purpose of storing intermediate information

		//	// If your protocol is to support periodic boundary conditions call this:
		//	getWSpace().cleanSpace(); // move any stray molecules back into simulation

		//	std::cout << "Step number:\t" << Step << std::endl;

		//	// Display your infoline every so often (UpdateScr)
		//	if((OutputLevel > Verbosity::Normal)&&(UpdateScr > 0)&&
		//		 (((Step) % UpdateScr) == 0))
		//	{
		//		infoLine();
		//	}
		//	

		//	// Save the coordinates in the trajectory as required
		//	if((OutputLevel > Verbosity::Normal)&&(UpdateTra > 0)&&
		//		 (((Step) % UpdateTra) == 0))
		//	{
		//		getWSpace().outtra.append();
		//	}
		//
		//}

		//// print some statistics if you want
		//if(OutputLevel > Verbosity::Normal) {
		//	printFinalStatistics();
		//}

		//// runcore should return the number of force/energy evaluations 
		//// done 
		//return Step;
	//}


	// This function prints a single line of information that is
	// printed every so often (every UpdateScr) steps to inform the user
	// about the state of the simulation.
	void BudeEMC::infoLine() const
	{
		// prints a line of current state of system, monitors and energies

		// Information on the state of the protocol, i.e. step number, temperature,
		// etc...
		printf("%5d",Step );  

		// Information on the current data in the monitors
		mon.printCurData();

		// Information on the current energies
		ff->infoLine();
		
		// Newline (the previous calls don't print a newline character)
		printf("\n");
	}

	// this function should "match up" and label the columns produced by the
	// above function infoLine() and is usually called at the beginning of a run
	void BudeEMC::infoLineHeader() const
	{
		// prints the headers for the above function
		printf("%5s", "Generation");

		// Headers of the monitors
		mon.printHeader();

		// Headers of the forcefield components
		ff->infoLineHeader();

		// New Line
		printf("\n");
	}

}
