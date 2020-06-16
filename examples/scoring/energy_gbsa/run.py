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

## and dy detail (useful for when things DO break)
ff.printEnergyByAtom()

## also show parameters
ff.info()


