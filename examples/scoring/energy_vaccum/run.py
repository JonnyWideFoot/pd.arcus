from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()
ffps.readLib("amber03aa.ff")

## test simple loading of PDB file
sim = PDB_In(ffps, "trpcage.pdb");
sim.loadAll();

# create workspace
wspace = WorkSpace( sim )
# print loaded system - this can be compared later
wspace.printPDB("inout.pdb")

# create a few common forcefields
ff = Forcefield(wspace)

# set up a simple vaccuum forcefield
bonds = BondedForcefield(wspace)
nb    = FF_NonBonded(wspace)
nb.Cutoff = 12.0
nb.InnerCutoff =  9.0
ff.add( bonds )
ff.add( nb )

## print energies as summary
ff.printEnergySummary()
## and dy detail (useful for when things DO break)
ff.printEnergyByAtom()
## also show parameters
ff.info()


