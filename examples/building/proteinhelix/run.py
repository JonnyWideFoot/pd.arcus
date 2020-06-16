from pd import *
cseed(4)
info()
timer()

ffps = FFParamSet()

## Use that Charmm forcefield
ffps.readLib( after("-ff","charmm22.ff")  )

## test simple loading of PDB file
sim = System( ffps );

## Create a Protein helix using evry amino acid
sim.add( NewProteinHelix(ffps,"*P-(ACDEFGHIKLMNPQRSTVWY)-R*") );


# create workspace
wspace = WorkSpace( sim )

# print loaded system - this can be compared later
wspace.printPDB("built.pdb")

# create a few common forcefields
ff = Forcefield(wspace)

# set up a simple vaccuum forcefielid and calculate the energy
bonds = FF_Bonded(wspace)
nb    = FF_NonBonded(wspace)
nb.Cutoff = 12.0
nb.InnerCutoff =  9.0
ff.add( bonds )
ff.add( nb )

## print energies as summary
ff.printEnergySummary()
ff.printEnergyByAtom()
## also show parameters
ff.info()


