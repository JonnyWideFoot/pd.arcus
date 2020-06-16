#include "global.h"

#include "../pd/accessory.h"

#include "system/fundamentals.h"

#include "system/genpolymer.h"
#include "workspace/workspace.h"

#include "sequence/sequence.h"
#include "sequence/alignment.h"

#include "forcefields/bude.h"
#include "forcefields/ffbonded.h"
#include "forcefields/nonbonded.h"
#include "forcefields/gbff.h"
#include "forcefields/lcpo.h"

#include "fileio/tra.h"
#include "fileio/trablocks.h"
#include "fileio/pdb.h"
#include "fileio/psfdcd.h"

#include "tools/io.h"
#include "tools/statclock.h"
#include "tools/stringtool.h"

#include "protocols/energy.h"
#include "protocols/minimise.h"
#include "protocols/torsionalminimisation.h"
#include "protocols/md.h"
#include "protocols/bude_emc.h"

#include "workspace/space.h"

#include "workspace/cluster.h"
#include "system/molecule.h"

using namespace std;
using namespace Physics;
using namespace Protocol;
using namespace Tra;
using namespace PDB;
using namespace Sequence;


void createSystem()
{
	FFParamSet ffps("../../pd/lib/bude01.ff");

	System sim(ffps);
	sim.add(NewProtein(ffps,"NALA-ALA-CALA"));
	sim.add(NewProtein(ffps,"NALA-ALA-CALA"));

	WorkSpace wspace( sim );
	wspace.printPDB("inpout.pdb");
}

void calcEnergies()
{
	FFParamSet ffps("../../pd/lib/bude01.ff");

	PDB_In dockingSystem(ffps,"4pep_prot.pdb");
	dockingSystem.setFilter(Library::Polypeptide); // ensure that we dont import waters and such
	dockingSystem.load(); // load all polypeptide from file 1 - the receptor

	PDB_In ligand(ffps,"pep4_A_84_237_mod.pdb");
	ligand.setFilter(Library::Polypeptide); // ensure that we dont import waters and such
	ligand.load(); // load all polypeptide from file 2 - the ligand

	dockingSystem.add(ligand); // add a copy of the ligand molecule to the 'dockingSystem'

	WorkSpace wspace( dockingSystem ); // create a WorkSpace

	Forcefield ff(wspace); // create a forcefield container

	BudeForcefield budeff(wspace); // create a BudeForcefield

	budeff.setCutoff(15.0);
	budeff.setCutoff_elec_formal(6.0);
	budeff.setCutoff_elec_partial(4.0);
	budeff.setReceptorIndex(0);
	budeff.setLigandIndex(1);

	ff.add( budeff ); // add the BudeForcefield to the forcefield container

	// print some info
	ff.printEnergySummary();
	ff.info();
	ff.infoLine();

	// set up printing to trajectories
	//Traj_Out tra("outtra",wspace);
	//wspace.addTra(tra);

	//PDB_Out pdbTra("pdbOutTra",wspace);
	//wspace.addTra(pdbTra);

	//OutTraPSF_DCD dcdOutput("dcd_output",&wspace);
	//wspace.addTra(dcdOutput);

}

void budeEMC()
{
	FFParamSet ffps("../../pd/lib/bude01.ff");

	//// You only need this once to convert the library format.
	//// I have included this for completeness.
	Library::RotamerLibrary rotLibTemp( ffps );
	rotLibTemp.convertLib("../../pd/lib/rotlib/bude.legacy.rotamer", Library::RotLibConvert_OldPDFormat() );
	rotLibTemp.writeLib( "../../pd/lib/rotlib/bude.rotamer" );

	// Every-day execution
	Library::RotamerLibrary rotLib( ffps );
	rotLib.readLib("../../pd/lib/rotlib/bude.rotamer" );

	PDB_In dockingSystem(ffps,"4pep_prot_chainA.pdb"); // load receptor
	dockingSystem.setFilter(Library::Polypeptide); // ensure that we dont import waters and such
	dockingSystem.load(); // load all polypeptide from file 1 - the receptor

	//Sequence::BioSequence seq(ffps);
	//seq.append("*C-(WC)-S*");

	PDB_In ligand(ffps,"pep4_A_84_237_ss.pdb"); // load ligand
	ligand.setFilter(Library::Polypeptide); // ensure that we dont import waters and such
	ligand.load(); // load all polypeptide from file 2 - the ligand
	//ligand.loadExplicit('I',Library::Polypeptide,seq); // load all polypeptide from file 2 - the ligand

	dockingSystem.add(ligand); // add a copy of the ligand molecule to the 'dockingSystem'

	WorkSpace wspace( dockingSystem ); // create a WorkSpace from the System

	wspace.printPDB("check_input.pdb");

	//return;

	Forcefield ff(wspace); // create a forcefield container

	BudeForcefield budeff(wspace); // create a BudeForcefield

	budeff.setCutoff(15.0);
	budeff.setCutoff_elec_formal(6.0);
	budeff.setCutoff_elec_partial(4.0);
	budeff.setReceptorIndex(0);
	budeff.setLigandIndex(1);

	ff.add( budeff ); // add the BudeForcefield to the forcefield container

	//// print some info
	//ff.printEnergySummary();
	//ff.info();
	//ff.infoLine();

	// set up printing to trajectories
	Traj_Out tra("outtra",wspace);
	wspace.addTra(tra);

	PDB_Out pdbTra("pdbOutTra",wspace);
	wspace.addTra(pdbTra);

	//OutTra_NAMD dcdOutput("dcd_output",&wspace);
	//wspace.addTra(dcdOutput);	

	// set up a budeEMC protocol
	BudeEMC budeEmc(ff,rotLib);

	budeEmc.setTransformations("trans_rot_777777.trans");
	budeEmc.setInitialGeneration("initial_pd_choices.txt");
	//budeEmc.setRotamerLibrary("../../pd/lib/bude_rotlib.pdb");
	budeEmc.setLigandRotamerPerturbations("rotchoice_ligand.rot");
	//budeEmc.setReceptorRotamerPerturbations("rotchoice_receptor.rot");
	//budeEmc.setRotamerChildrenCount(105);

	budeEmc.setInitialRotamerChildrenCount(10);
	budeEmc.setRotamerGenerationSize(12,5); // 12 Parents, 5 Children per parent
	budeEmc.setRotamerSteps(3);

	budeEmc.RandomSeed = 86121689;
	budeEmc.InitialGenCount = 20;
	budeEmc.ChildrenPerParent = 7; // 7 children per parent for 15 parents gives generation size of 105.
	budeEmc.MutationRate = 0.8;
	budeEmc.NextGenCount = 15;

	budeEmc.ReceptorIndex = 0;
	budeEmc.LigandIndex = 1;
	budeEmc.Steps = 3; // == number of generations
	budeEmc.UpdateScr = 1;
	budeEmc.OutputLevel = Verbosity::Normal;
	budeEmc.run();

}

void cluster()
{
	FFParamSet ffps("../../pd/lib/bude01.ff");

	System sim(ffps);
	sim.add(NewProtein(ffps,"NALA-ALA-CALA"));
	sim.add(NewProtein(ffps,"NALA-ALA-CALA"));

	WorkSpace wspace( sim );
	wspace.printPDB("inpout.pdb");

	Forcefield ff(wspace);
	BondedForcefield bonds(wspace);
	FF_NonBonded nb(wspace);
	nb.Cutoff =  9.0;
	nb.InnerCutoff =  6.0;
	ff.add( bonds );
	ff.add( nb );

	Minimisation min(ff);
	min.Steps = 100;
	min.UpdateScr  = 1 ;
	min.UpdateTra  = 0;
	min.run();

	MolecularDynamics md(ff);
	md.Steps = 1000;
	md.Timestep   = 1.00E-15;
	md.UpdateScr = 100;
	md.UpdateTra = 100;
	md.UpdateNList = 10;
	md.TargetTemp = new Temp(3000);


	ClusterSimple myCluster1;

	SnapShotOutTra snapout(myCluster1,wspace);
	wspace.addTra(snapout);

	Traj_Out clustered("clustered",wspace);
	TraBlockCmnt comment(256);
	clustered.addBlock(comment);


	myCluster1.clusterCutoff = 2;

	for(int i = 0; i < 5; i++)
	{
		md.run(); // run another set of steps
		//myCluster1.push_back( psystem.save() );
		//snapout.append();
	}

	myCluster1.doClustering();

	std::cout << "Number of index groups:\t" << myCluster1.clusterSize() << std::endl;

	for(int i = 0; i < myCluster1.clusterSize(); i++)
	{
		std::cout << "Center of cluster number\t " << i << ":\t" << myCluster1.getCentre(i) << std::endl;

		std::cout << "Contents of cluster number\t " << i << ":\t";

		for(int j = 0; j < myCluster1.indexSize(i); j++)
		{			
			const SnapShot &data = myCluster1.getData(myCluster1.getIndexMember(i,j));
			wspace.load(data);

			std::string report;
			report = "Contents of cluster number\t " + int2str(i) + ":\t" + int2str(myCluster1.getIndexMember(i,j));
			comment.setTo(report.c_str());
			clustered.append();

			std::cout << myCluster1.getIndexMember(i,j) << ",";

		}
	}
}

int main()
{
	//createSystem()

	//calcEnergies();
	
	budeEMC();

	//cluster();

	return 0;

}



