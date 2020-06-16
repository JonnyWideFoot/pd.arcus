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
bonds = BondedForcefield(wspace)
gb    = GB_Still( wspace)
gb.Cutoff = 12.0
gb.InnerCutoff =  9.0
gb.FastMode = 1
sasa  = SASA_LCPO( wspace )
sasa.GlobalASP = 0.005

## add all the forcefield components to our forcefield
ff.add( bonds )
ff.add( gb )
ff.add( sasa )

## print energies as summary
ff.printEnergySummary()

## also show parameters
ff.info()

## try some minimisation 
min = Minimisation(ff)
min.Steps = 25
min.UpdateScr  = 1 
min.UpdateTra  = 0
min.run()

## restore and try an MD run (constant TEMP)
## this should keep its temperature 
md = MolecularDynamics(ff)
md.Steps = 50 
md.Timestep   = 1.00E-15
md.UpdateScr = 10
md.UpdateTra = 100
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.LangevinOnHydrogens = True
md.setTemperature(300)


rex = REX_Local(ff);
rex.Steps = 4;
rex.OutputLevel = Verbosity.Loud
rex.addReplica( md, 300 )
rex.addReplica( md,  320 )

rex.run()

timer()

