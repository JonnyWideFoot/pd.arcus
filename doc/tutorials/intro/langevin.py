from pd import *
ctimeseed()   ## set random number seed to current system time
timer()       ## start a timer 


ffps = FFParamSet("../../../lib/amber03aa.ff")

## Create a simulation system using the
## parameter set ffps
sim = System(ffps)
sim.add(NewProteinPDB(ffps,"trpcage.pdb"))
sim.info()

## Create a simulation workspace
wspace = WorkSpace( sim )
wspace.info()

## Create an output trajectory (tra format) 
tra = OutTraPSF_DCD("output",wspace)
wspace.addTra(tra)

## create a GB/SA Forcefield using the AMBER
## 03 Parameter set

ff = Forcefield(wspace)

bonds = BondedForcefield()

gb    = GB_Still()
gb.Cutoff = 12.0         ## Cutoffs in Angstrom
gb.InnerCutoff =  9.0
gb.VdwCutoff =  8.0
gb.VdwInnerCutoff =  6.0
gb.FastMode = False #True 

sasa  = SASA_LCPO()
sasa.GlobalASP = 0.005   ## surface tension +5cal/mol

## Add forcefield components to forcefield 
ff.add( bonds )
ff.add( gb )
ff.add( sasa )


ff.printEnergySummary()  ## calculate energies and show them summarised
ff.info()         ## show parameters


## Run some Minimisation
min = Minimisation(ff)
min.Steps = 50
min.UpdateScr  =  10
min.UpdateNList = 10
min.run()

## Run some Molecular Dynamics
md = MolecularDynamics(ff)
md.Steps = 2000
md.UpdateScr = 100
md.UpdateTra = 10
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.TargetTemp = Temp(300)

# Add a atom-atom distance monitor to the simulation 
# (we will monitor atoms with indices 102 and 263)
dist = AtomDistanceMonitor(wspace, 102, 263)
md.addMonitor(dist)
md.run()

## put the collected data into a histogram with bin width 0.05
dist.getHistogram(0.05).print_count()

## print time taken
timer()

