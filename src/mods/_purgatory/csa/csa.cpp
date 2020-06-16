// --------------------------------------------------------------------------------
// File: csa.cpp
// Description: Conformational Space Annealing (Scheraga et al. 1999)
// Version: 0.1
// Author: Michael Tyka
// mail: m.tyka@bristol.ac.uk
//
// Dependencies: wspace.h
// Headers: csa.h
// Sources: csa.cpp
// Classe(s): ConformationalSpaceAnnealing

// Pre-compiled Header
#include "global.h"

// Self Header Include
#include "csa.h"

// mmlib includes
#include "protocols/minimise.h"
#include "protocols/dualffminimiser.h"
#include "protocols/md.h"
#include "protocols/montecarlo.h"
#include "fileio/tra.h"
#include "library/angleset.h"
#include "workspace/workspace.h"

using namespace Maths;
using namespace Physics;
using namespace Protocol;

namespace Protocol
{
	class IndexBasedGroup{
	public:
		IndexBasedGroup():member(200) {
			lowestEnergyIndex = -1;
		};
		~IndexBasedGroup(){
		}

		IndexBasedGroup &operator=(const IndexBasedGroup & ibg){
			member = ibg.member;
			lowestEnergyIndex = ibg.lowestEnergyIndex;
			lowestEnergy = ibg.lowestEnergy;
			return (*this);
		}

		std::vector < int >member;
		int lowestEnergyIndex;
		double lowestEnergy;

		int findExtreme_Lowest(double *property){
			double extProp = (double)(1E30);
			int extI = -1;
			for(int i = 0; i < member.size(); i++) {
				if(property[member[i]] < extProp) {
					extProp = property[member[i]];
					extI = i;
				}
			}
			return extI;
		}
		int findExtreme_Highest(double *property){
			double extProp = (double)(-1E30);
			int extI = -1;
			for(int i = 0; i < member.size(); i++) {
				if(property[member[i]] > extProp) {
					extProp = property[member[i]];
					extI = i;
				}
			}
			return extI;
		}

		void sortBy(double *property){
			int *indexarray = new int[member.size()];
			double *proparray = new double[member.size()];
			int i;

			for(i = 0; i < member.size(); i++)
				indexarray[i] = member[i];
			memcpy(proparray, property, sizeof(double) * member.size());

			qcksort(proparray, indexarray, member.size());

			// write indices back into dynamic array
			for(i = 0; i < member.size(); i++)
				member[i] = indexarray[i];

			delete[]indexarray;
			delete[]proparray;
		}
	};


	class PSPIterativeGrouping{
	public:
		PSPIterativeGrouping(WorkSpace * _wspace,
			SnapShot * _psparray, int _nstructures, double RMScutoff, int _maxgroupsize){
				energylist = NULL;
				sortedEnergyList = NULL;
				sortedEnergyListIndex = NULL;
				nstructures = _nstructures;
				psparray = _psparray;
				wspace = _wspace;
				crmsgrouplimit = RMScutoff;
				maxgroupsize = _maxgroupsize;
		}

		friend class ConformationalSpaceAnnealing;

		~PSPIterativeGrouping(){
			delete[]energylist;
			delete[]sortedEnergyList;
			delete[]sortedEnergyListIndex;
		}

		int prepareGrouping(){
			int i;
			group.clear();

			delete[]energylist;
			delete[]sortedEnergyList;
			delete[]sortedEnergyListIndex;

			energylist = new double[nstructures];
			sortedEnergyList = new double[nstructures];
			sortedEnergyListIndex = new int[nstructures];

			printf("Gathering potential energies .. \n");
			for(i = 0; i < nstructures; i++) {
				energylist[i] = psparray[i].epot;
				sortedEnergyList[i] = psparray[i].epot;
				sortedEnergyListIndex[i] = i;
			}

			printf("Sorting potential energies .. \n");
			qcksort(&sortedEnergyList[0], &sortedEnergyListIndex[0], nstructures);
			return 0;
		}

		int initialGrouping(double cRMSlimit){
			// normal variables
			int n, i;

			/*
			// for(i=0;i<100;i++){
			int j=1;
			double j_ene = 20;
			IndexBasedGroup newgroup;
			newgroup.member.push_back(&j);
			newgroup.lowestEnergy = j_ene;
			newgroup.lowestEnergyIndex = j;
			group.push_back(&newgroup);
			// }

			IndexBasedGroup *ibg = new IndexBasedGroup[100];

			for(i=0;i<100;i++) ibg[i] = newgroup;

			delete [] ibg;

			return 0;
			*/
			printf("Generating groups .. (from %5d structures) \n", nstructures);
			// now start adding structures one by one in order of increasing energy
			for(i = 0; i < nstructures; i++) {
				int lowestCluster = -1;
				double lowestcRMS = 100000;
				int j = sortedEnergyListIndex[i]; // current structure
				double j_ene = sortedEnergyList[i]; // Current energy

				if((j_ene > 1500) || // dont include structures with outrageous energies (probably corrupted)
					(j_ene < -1500))
					continue;

				int k; // group representative
				for(n = 0; n < group.size(); n++) {
					k = group[n].lowestEnergyIndex; // get group representative index

					double cRMS = psparray[j].cRMSFrom(psparray[k], *wspace, "CA");
					if(cRMS < lowestcRMS) {
						lowestCluster = n;
						lowestcRMS = cRMS;
					}
				}

				if((lowestCluster < 0) || (lowestcRMS > cRMSlimit)) {// new cluster
					// printf("Creating new cluster %d for %d \n",group.size()+1,j);
					IndexBasedGroup newgroup;
					newgroup.member.push_back(j);
					newgroup.lowestEnergy = j_ene;
					newgroup.lowestEnergyIndex = j;
					group.push_back(newgroup);
				} else {// add to existing cluster
					//printf("adding %d to existing cluster \n",j);


					group[lowestCluster].member.push_back(j);// add index j to this group
					// if new energy is lower, change group representative (shouldn't happen though)
					// since everything is done in order of increasing energy anyway...
					if(group[lowestCluster].lowestEnergy > j_ene) {
						printf("Strange ? \n");
						group[lowestCluster].lowestEnergy = j_ene;
						group[lowestCluster].lowestEnergyIndex = j;
					}
				}

				//delete [] jtraCA;
			}




			return 0;
		}

		/// --------------------------------------------------------------------

		int refineGrouping(double cRMSlimit){
			/* int i,n;
			printf("Refining groups\n");
			int refchanges=0;
			int igroup;
			int imember;

			for(igroup=0;igroup<group.size();igroup++){
			for(imember=0;imember<group[igroup]->member.size();imember++){
			int lowestCluster = -1;
			double lowestcRMS = 100000;
			int j = *(group[igroup]->member[imember]); // current structure
			double j_ene = energylist[j]; // Current energy
			int k; // group representative

			double owncrms;
			//printf("%6d (%6d) --> ",i,j);
			for(n=0;n<group.size();n++){
			k = group[n]->lowestEnergyIndex; // get group representative index
			double cRMS = calccRMS(&alltraCA[nAtoms*j],&alltraCA[nAtoms*k],nAtoms);
			if(n==igroup) owncrms = cRMS;
			if(cRMS < lowestcRMS){
			lowestCluster = n;
			lowestcRMS = cRMS;
			}
			}

			if(lowestCluster == igroup){//are we still best in our own group ?
			// printf("No change required\n");
			continue;
			}else{
			refchanges++;
			//printf("Imp: %lf %lf\n",owncrms,lowestcRMS);
			// first remove from old group

			if(*group[igroup]->member[imember] ==
			group[igroup]->lowestEnergyIndex){
			printf("Removing centre ?? \n");
			}
			group[igroup]->member.remove(imember);

			// and add to best new group (if
			if((lowestCluster < 0)||(lowestcRMS > cRMSlimit) ){ // new cluster
			printf("Creating new cluster %d for %d \n",group.size()+1,j);
			GenericIndexBasedGroup newgroup;
			newgroup.member.push_back(&j);
			newgroup.lowestEnergy = j_ene;
			newgroup.lowestEnergyIndex = j;
			group.push_back(&newgroup);
			}else{ // add to existing cluster
			// printf("adding %d to existing cluster \n",j);
			group[lowestCluster]->member.push_back(&j); // add index j to this group
			// if new energy is lower, change group representative (shouldn't happen though)
			// since everything is done in order of increasing energy anyway...
			if(group[lowestCluster]->lowestEnergy > j_ene){
			group[lowestCluster]->lowestEnergy = j_ene;
			group[lowestCluster]->lowestEnergyIndex = j;
			}
			}
			}

			}
			}
			*/
			// printf("Refinement changes : %d\n",refchanges);
			// return refchanges;
			return 0;
		}


		/// --------------------------------------------------------------------

		int runcore(){
			int starttime = (int) time(NULL);
			if(prepareGrouping() != 0)
				return -1;
			if(initialGrouping(crmsgrouplimit) != 0)
				return -1;
			// for(int refinement=0;refinement < 50; refinement++){
			// if(refineGrouping(crmsgrouplimit)<=0) break;
			// }
			int endtime = (int) time(NULL);
			printf("Clustering time: %d seconds \n", endtime - starttime);
			return 0;
		}

	private:
		std::vector < IndexBasedGroup > group;

		double *energylist;
		double *sortedEnergyList;
		int *sortedEnergyListIndex;

		double crmsgrouplimit;
		int maxgroupsize;

		SnapShot *psparray;
		int nstructures;

		WorkSpace *wspace;

	};







	void ConformationalSpaceAnnealingPS::info(){
		printf("csa.param.Steps %d\n", Steps);
		printf("csa.param.librarysize %d\n", librarysize);
		printf("csa.param.maxgroupsize %d\n", maxgroupsize);
		printf("csa.param.Temperature %lf\n", Temperature);
		printf("csa.param.mcmsteps %d\n", mcmsteps);
		printf("csa.param.mcmTemperature %lf\n", mcmTemperature);
		printf("csa.param.mcmpreminsteps %d\n", mcmpreminsteps);
		printf("csa.param.mcmminsteps %d\n", mcmminsteps);
		printf("csa.param.rmslimit %lf\n", rmslimit);
		printf("csa.param.rmslimit_end %lf\n", rmslimit_end);
		printf("csa.param.exchange_prop: %lf\n", exchange_prop);
	}



	//--------------------------------------------------------------
	// Begin: Conformational Space Annealing
	//--------------------------------------------------------------

	ConformationalSpaceAnnealing::ConformationalSpaceAnnealing(WorkSpace * wspace,
		Physics::Forcefield * _ff,
		Physics::Forcefield * _stericff,
		Library::AngleSet * _angset,
		std::vector < Manipulator::MoveBase * >* _sp_csa,
		std::vector < Manipulator::MoveBase * >* _sp_fine,
		ConformationalSpaceAnnealingPS &_ps):
	ProtocolBase(wspace,_ff),
		sp_csa(_sp_csa),
		sp_fine(_sp_fine),
		stericff(_stericff),
		angset(_angset),
		ps(_ps)
	{
		psplib = new SnapShot[ps.librarysize];
	};

	ConformationalSpaceAnnealing::~ConformationalSpaceAnnealing(){
		delete[] psplib;
	};

	void ConformationalSpaceAnnealing::info() const{ // prints a little block of parameter information
	}
	void ConformationalSpaceAnnealing::infoLine() const { // prints a line of current energies
		int test = Particle::flag_Backbone;
	}
	void ConformationalSpaceAnnealing::infoLineHeader() const { // prints the headers for the above function
	}

	int ConformationalSpaceAnnealing::runcore()
	{
		int ilib;
		int ipert, pertdone;
		int move;
		int nminimisations = 0;
		double lowestEnergy = 100000;
		double lowestCRMS = 100;
		char statfilename[100];
		char flibfilename[100];
		char flib2ilename[100];

		int i,j;

		sprintf(&statfilename[0], "%s.csa.stat", ps.stem.c_str());
		i = rand() % 100;
		sprintf(&flibfilename[0], "%s.csa.tra", ps.stem.c_str());
		sprintf(&flib2ilename[0], "%s.csasa.tra", ps.stem.c_str());

		int newgen_memorysize = ps.librarysize * (ps.maxgroupsize + ps.repeats + 1);
		SnapShot *newgeneration = new SnapShot[newgen_memorysize];
		SnapShot *survivor = new SnapShot[newgen_memorysize];
		SnapShot *finallib = new SnapShot[ps.librarysize];
		int finallibsize = 0;
		int nnewgen;

		printf("CSA Parameters: \n");
		printf("stat filename: %s\n", &statfilename[0]);
		printf("csa.tra filename: %s\n", &flibfilename[0]);
		ps.info();
		printf("sizeof PSP: %d bytes\n", sizeof(SnapShot) +
			sizeof(SnapShotAtom) *
			wspace->atom.size());
		printf("Small pertubations: %d \n", sp_fine->size());
		printf("Large pertubations: %d \n", sp_csa->size());

		// initialise all the PSPs
		for(ilib = 0; ilib < newgen_memorysize; ilib++) {
			newgeneration[ilib] = wspace->save();
		}
		nnewgen = 0;

		// add the current starting structures (stored in psplib) to the new generation
		for(ilib = 0; ilib < ps.librarysize; ilib++) {
			newgeneration[nnewgen] = wspace->save();
			newgeneration[nnewgen] = psplib[ilib];
			nnewgen++;
		}

		for(Step = 0; Step < ps.Steps; Step++) {
			printf("CSA Step: %d \n", Step);
			double library_lowestCRMS = 10000;

			for(ilib = 0; ilib < ps.librarysize; ilib++) {

				wspace->load(psplib[ilib]);

				if(frand() < ps.exchange_prop) {
					int startir;
					int endir;
					int istruc = rand() % ps.librarysize;
					SnapShot pspx;
					pspx = wspace->save();

					for(int htry = 0; htry < 10; htry++) {
						startir = rand() % (wspace->res.size() - 4);
						endir = startir + rand() % 7 + 3;
						if(endir > (wspace->res.size() - 1))
							endir = wspace->res.size() - 1;
						if(pspx.insertPSPfragment(psplib[istruc], startir, endir, *wspace) == 0)
							break;
					}
					wspace->load(pspx);
					//printf("inserted fragment %d<-->%d from %d \n",startir,endir,istruc);
					pertdone = -1;
				} else {
					pertdone = 0;

					//for(ipert=0;ipert<sp_csa->size();ipert++){

					// large move
					if(sp_csa->size() > 0) {
						do {
							ipert = rand() % sp_csa->size();
							move = ((*sp_csa)[ipert])->apply();

							if(move != 0)
								pertdone = ipert + 1; // note down which pertubation was done
						} while(pertdone == 0);
					}
					//small move
					for(ipert = 0; ipert < sp_fine->size(); ipert++) {
						move = ((*sp_fine)[ipert])->apply();
					}
				}
				//wspace->outtra.append();

				if((ps.mcmsteps == 0) && (ps.mcmminsteps == 0) && (ps.mcmpreminsteps == 0)) {
					nminimisations += 1;
					wspace->nlist().calcNewList();
					ff->calcEnergies();
				} else {
					DualFFMinimiser dualmin(wspace,ff, ps.mcmminsteps, stericff, ps.mcmpreminsteps, 0.05 * PhysicsConst::kcal2J / PhysicsConst::Na);
					Temp mcTemp(ps.mcmTemperature);
					MonteCarlo mc(wspace,ff,&dualmin,ps.mcmsteps,&mcTemp,sp_fine);
					nminimisations += ps.mcmsteps;
					mc.UpdateTra = 0;
					mc.Silent = true;
					if(mc.runcore() < 0) // if it fails, i.e. returns -1, not a number of Steps
						continue;
					wspace->load(mc.lowestEnergyPSP);
				}

				wspace->ene.data1 = (double)pertdone;

				// wspace->stdff->infoLine();
				// wspace->outtra.append();

				int ngentry;
				if(Step > 0) { // unless first round, first repeat add to library
					ngentry = nnewgen;
					nnewgen++;
				} else { // otherwise replace library
					ngentry = ilib;
				}

				if(ngentry >= newgen_memorysize) { // this should in theory not happen (newgen_memorysize should be sufficiently large) but just in case to prevent crashes check anyway
					printf("CODE ERROR: newgenration[] overflow (Step=%d,ilib=%d) - continueing \n", Step, ilib);
					ilib = ps.librarysize;
					break;
				}
				newgeneration[ngentry] = wspace->save();
				newgeneration[ngentry].status2 = Step;// record it's creation time;
				// newgeneration[ngentry].setTo(mcm.lowestEnergyPSP,false);
				printf("--> %4d %8.3lf %d\n", ilib, double (newgeneration[ngentry].epot) * PhysicsConst::J2kcal * PhysicsConst::Na, pertdone);

				fflush(stdout);
			}

			// double cRMS = newgeneration[0].dRMSFrom(&newgeneration[1],0,wspace);
			double rmslimit_step = (ps.rmslimit_end * (double) Step / (double) ps.Steps) +
				(ps.rmslimit * (1.0 - (double) Step / (double) ps.Steps));
			printf("Grouping procedure: RMS limit: %lf \n", rmslimit_step);
			PSPIterativeGrouping pspig(wspace, &newgeneration[0], nnewgen, rmslimit_step, 4);
			pspig.runcore();

			printf("Refilling Library (size=%d)\n", ps.librarysize);

			fflush(stdout);

			//refill library
			int igroup = 0;
			int imember;
			int ichoice;

			for(ilib = 0; ilib < ps.librarysize; ilib++) {
				double *problist = new double[pspig.group[igroup].member.size()];

				//printf("Ilib %d: Choosing probability \n",ilib);
				for(imember = 0; imember < pspig.group[igroup].member.size(); imember++) {
					problist[imember] = exp(-(double (newgeneration[pspig.group[igroup].member[imember]].epot) -
						double (newgeneration[pspig.group[igroup].lowestEnergyIndex].epot)) /(PhysicsConst::kB *
						ps.Temperature));
				}

				//choose next candidate by boltzmann distributution
				ichoice = chooseByProbability(&problist[0], pspig.group[igroup].member.size());
				//ichoice = pspig.group[igroup].lowestEnergyIndex;

				//load this choice
				psplib[ilib] = newgeneration[pspig.group[igroup].member[ichoice]];

				igroup++;
				if(igroup >= pspig.group.size())
					igroup = 0; // if we have arrived at the last group, start from the beginning

				delete[]problist;
			}

			fflush(stdout);
			//display groups
			for(i=0;i<pspig.group.size();i++){
				if(i<ps.librarysize)printf("+ "); else printf("- ");
				printf("%4d %8.2lf\n",pspig.group[i].lowestEnergyIndex,
					double(newgeneration[pspig.group[i].lowestEnergyIndex].epot) 
					* PhysicsConst::J2kcal * PhysicsConst::Na);
			}

			printf("Refilling survivors library...");
			// refill generation bank
			int nsurvivors = 0;

			for(i = 0; i < min((size_t)ps.librarysize, pspig.group.size()); i++) {
				// the best 'maxgroupsize' in each group survive

				// pspig.group[i]->sortBy(pspig.energylist);
				// not needed because it's already sorted (by the clustering algorithm)

				int leaderpersistance;
				for(int m = 0; m < min((size_t)ps.maxgroupsize, pspig.group[i].member.size()); m++) {
					int structureIndex = pspig.group[i].member[m];
					survivor[nsurvivors] = wspace->save();
					survivor[nsurvivors] = newgeneration[structureIndex];
					if(m == 0) { // if it's the group leader then add up its leader number (status3)
						survivor[nsurvivors].status3++;
						leaderpersistance = survivor[nsurvivors].status3;

						wspace->load(survivor[nsurvivors]);
						wspace->calccRMS();

						if(wspace->ene.cRMS < lowestCRMS)
							lowestCRMS = wspace->ene.cRMS;
						if(wspace->ene.cRMS < library_lowestCRMS)
							library_lowestCRMS = wspace->ene.cRMS;

						if(survivor[nsurvivors].status4 == 0) {
							survivor[nsurvivors].status4 = 1; // mark that this structure has been saved
							wspace->Step = Step;
							wspace->ene.data1 = (double)survivor[nsurvivors].status2; // creation Step
							wspace->ene.data2 = (double)survivor[nsurvivors].status3; // group leader age
							if(ps.UpdateTra > 1)
								wspace->outtra.append();
						}

						if(Step >= (ps.Steps - 1)) { // at the end of the run
							finallib[finallibsize] = wspace->save();
							finallibsize++;
						}

					} else {
						survivor[nsurvivors].status3 = leaderpersistance;
					}
					nsurvivors++;
				}
			}

			// copy over to new generation list
			for(i = 0; i < nsurvivors; i++) {
				newgeneration[i] = survivor[i];
			}
			nnewgen = nsurvivors;

			if(double (psplib[0].epot) < lowestEnergy)
				lowestEnergy = double (psplib[0].epot);

			printf("Newgen: %d \n", nnewgen);
			FILE *file = fopen(&statfilename[0], "a");
			fprintf(file, "%10d %9.3lf %9.3lf %9.3lf\n", nminimisations, lowestEnergy * PhysicsConst::J2kcal * PhysicsConst::Na, lowestCRMS,
				library_lowestCRMS);
			fclose(file);
		}


		// now analyse the resulting population ..
		printf("Analysing final library contents\n");
		printf("Filename: %s \n", &flibfilename[0]);
		double final_lowestcRMS = 1000000;
		int final_lowestcRMS_index = -1;
		double final_lowestdRMS = 1000000;
		double final_avene = 0;

		if(ps.UpdateTra > 0) {

			std::string filename(&flibfilename[0],strlen(&flibfilename[0]));
			Tra::Traj_Out ftra(filename,wspace);
			ftra.create(Tra::Energies);

			for(i = 0; i < finallibsize; i++) {
				wspace->load(finallib[i]);
				final_avene += finallib[i].epot;
				wspace->calccRMS();

				printf("csafinal %4d\t%6.2lf\t%8.3lf\n", i, wspace->ene.cRMS, double (wspace->ene.epot) * PhysicsConst::J2kcal * PhysicsConst::Na);

				if(wspace->ene.cRMS < final_lowestcRMS) {
					final_lowestcRMS = wspace->ene.cRMS;
					final_lowestcRMS_index = i;
				}
				if(wspace->ene.dRMS < final_lowestdRMS)
					final_lowestdRMS = wspace->ene.dRMS;

				ftra.append();
			}
		}
		final_avene /= finallibsize;
		printf("csafinal.lowene %8.2lf \n", lowestEnergy * PhysicsConst::J2kcal * PhysicsConst::Na);
		printf("csafinal.lowcrms %6.2lf %8.2lf %4d %5.3lf\n", final_lowestcRMS,
			double (finallib[final_lowestcRMS_index].epot) * PhysicsConst::J2kcal * PhysicsConst::Na,
			final_lowestcRMS_index, (double) final_lowestcRMS_index / (double) ps.librarysize);
		printf("csafinal.lowdrms %6.2lf \n", final_lowestdRMS);
		printf("csafinal.avene %6.2lf \n", final_avene * PhysicsConst::J2kcal * PhysicsConst::Na);

		wspace->load(finallib[0]); // load lowest lowest lowest conformation


		if(ps.finalsasteps > 0){
			printf("Simulated annealing of final library... \n");
			for(i = 0; i < finallibsize; i++) {
				printf("finallib %5d\n",i);

				wspace->load(finallib[i]);
				// sort out Temperature control behaviour (only has an effect if a Thermostat is switched on at the same time)

				MolecularDynamics mymd(wspace,ff);
				mymd.Integrator = MolecularDynamics::Langevin;
				mymd.Timestep = 2E-15;
				mymd.Steps = ps.finalsasteps;

				mymd.TargetTemp = new TempLinear(ps.finalsastartt,ps.finalsaendt);
				mymd.Thermostat = MolecularDynamics::NoThermostat;
				mymd.FricCoeff = (10.0*1E12);
				mymd.UpdateScr = 1000;
				mymd.UpdateTra = 0;
				mymd.Silent = false;
				mymd.runcore();

				delete mymd.TargetTemp;

				if(ps.finalsaminsteps>0){
					Minimisation minimise(wspace, ff);
					minimise.Algorithm = Minimisation::ConjugateGradients;
					minimise.Steps = ps.finalsaminsteps;
					minimise.StepSize = 0.1E7;
					minimise.UpdateScr = 100;
					minimise.Silent = false;
					minimise.runcore();
				}
				finallib[i] = wspace->save();
			}


			// now save the new structures
			if(ps.UpdateTra > 0) {
				final_lowestcRMS = 1000000;
				final_lowestcRMS_index = -1;
				final_lowestdRMS = 1000000;
				lowestEnergy = DBL_MAX;
				final_avene = 0;

				std::string filename(&flib2ilename[0],strlen(&flib2ilename[0]));
				Tra::Traj_Out ftra(filename,wspace);
				ftra.create(Tra::Energies);

				// sort by ene & print & put into a tra file

				double *finalsa_ene = new double[finallibsize];
				int *finalsa_ene_index = new int[finallibsize];

				for(i = 0; i < finallibsize; i++) {
					finalsa_ene[i] = finallib[i].epot;
					finalsa_ene_index[i] = i;
				}

				qcksort(finalsa_ene,finalsa_ene_index,finallibsize);
				for(j = 0; j < finallibsize; j++) {
					i = finalsa_ene_index[j];

					wspace->load(finallib[i]);
					final_avene += finallib[i].epot;
					if(finallib[i].epot < lowestEnergy){
						lowestEnergy = finallib[i].epot;
					}

					wspace->calccRMS();

					printf("csafinalsa %4d\t%6.2lf\t%8.3lf\n", i, wspace->ene.cRMS, double (wspace->ene.epot) * PhysicsConst::J2kcal * PhysicsConst::Na);

					if(wspace->ene.cRMS < final_lowestcRMS) {
						final_lowestcRMS = wspace->ene.cRMS;
						final_lowestcRMS_index = i;
					}
					if(wspace->ene.dRMS < final_lowestdRMS){
						final_lowestdRMS = wspace->ene.dRMS;
					}
					ftra.append();
				}
			}

			final_avene /= finallibsize;
			printf("csafinalsa.lowene %8.2lf \n", lowestEnergy * PhysicsConst::J2kcal * PhysicsConst::Na);
			printf("csafinalsa.lowcrms %6.2lf %8.2lf %4d %5.3lf\n", final_lowestcRMS,
				double (finallib[final_lowestcRMS_index].epot) * PhysicsConst::J2kcal * PhysicsConst::Na,
				final_lowestcRMS_index, (double) final_lowestcRMS_index / (double) ps.librarysize);
			printf("csafinalsa.lowdrms %6.2lf \n", final_lowestdRMS);
			printf("csafinalsa.avene %6.2lf \n", final_avene * PhysicsConst::J2kcal * PhysicsConst::Na);

			wspace->load(finallib[0]); // load lowest lowest lowest conformation


		} // if

		delete[]newgeneration;
		delete[]survivor;
		delete[]finallib;

		return Step;
	}

	int ConformationalSpaceAnnealing::fillLibrary(Library::RotamerLibrary *rotlib){
		int i;
		BackbonePhiPsiSetMove bppsp(wspace, rotlib, angset, 0.0, 1, false, false);

		for(i = 0; i < ps.librarysize; i++) {
			if(ps.scramble != 0)
				bppsp.changeConsResiduesRandomly(0, wspace->res.size());
			psplib[i] = wspace->save();
		}
		return 0;
	}
} // namespace 'Protocol'


