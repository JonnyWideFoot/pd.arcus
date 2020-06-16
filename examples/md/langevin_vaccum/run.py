from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()
ffps.readLib("amber03aa.ff")

## test simple loading of PDB file
sim = PDB_In(ffps, "trpcage.pdb");
sim.loadAll();

# create workspace from the system sim
wspace = WorkSpace( sim )

# print loaded system - this can be compared later
wspace.printPDB("inout.pdb")

# create a few common forcefield components

ff = Forcefield(wspace)
# set up a simple vaccuum forcefield
bonds = BondedForcefield(wspace)
nb    = FF_NonBonded(wspace)
nb.Cutoff = 20.0
nb.InnerCutoff =  16.0
ff.add( bonds )
ff.add( nb )

## add all the forcefield components to our forcefield
ff.add( bonds )
ff.add( gb )
ff.add( sasa )

## print energies as summary
ff.printEnergySummary()

## also show parameters
ff.info()


## preminimise to prevent the MD to explode straight away 
min = Minimisation(ff)
min.Steps = 25
min.UpdateScr  = 1 
min.UpdateTra  = 0
min.run()


## Set up Langevin Dynamics which is handles by the class MolecularDynamics

md = MolecularDynamics(ff)   ## create the MD protocol
md.Steps = 2000              ## set the number of simulation steps
md.Timestep   = 0.5E-15      ## set the timestep to 1 femto second             
md.UpdateScr = 10            ## Print a line to the screen every 10 steps
md.UpdateTra = 50            ## Save structures to each trajectory every 50 steps
md.UpdateNList = 10          ## Update the Neighbour List every 10 steps

## Specific Langevin related settings
## NOTE: No thermostat should be set here ! Langevin dynamics has temperature control
## built in so to speak, due to the way it balances random and friction forces.

## Set Langevin dynamics to be ON
md.Integrator = MolecularDynamics.Langevin

md.FricCoeff = 10E12         ## set friction coefficient for langevin 
                             ## dynamics (to 10E12 per second or 10 per picosecnd)							
md.LangevinOnHydrogens = True ## Do we want Hydrogens to be treated with stochastic 
                             ## simulation or using velcoity verlet (newtonian)
md.TargetTemp = Temp(300)    ## Run Langevin at 300K 
md.run()

ff.printEnergySummary()

timer()

