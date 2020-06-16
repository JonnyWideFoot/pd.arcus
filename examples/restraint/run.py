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

# positionally restrain all atoms of residue 0
crest = FF_Restraint_Positional(wspace);
crest.setSelection( PickResidue( 0 ) ); 
crest.k = 20.0   ## restraint constant in kcal/mol/A^2
crest.detail()
ff.add(crest);

# positionally restrain all atoms of residue 0
irest = FF_Restraint_Internal(wspace);
irest.setSelection( PickMolecule( 0 ) ); 
irest.k = 20.0   ## restraint constant in kcal/mol/A^2
irest.detail()
ff.add(irest);

# positionally restrain all atoms of residue 0
trest = FF_Restraint_Torsional(wspace);
trest.OneRestraintPerBond = True
trest.setSelection( PickResidue( 0 ) ); 
trest.k = 0.01   ## restraint constant in kcal/mol/rad
trest.info()
trest.detail()
ff.add(trest);

# special native contact restraint 
ncrest = FF_Restraint_NativeContact(wspace);
ncrest.setSelection( PickBackbone() ); 
ncrest.k = 1000   ## restraint constant in kcal/mol/rad
ncrest.info()
ncrest.detail()
ff.add(ncrest);

# special native contact restraint 
atrest = FF_Restraint_AtomDistance(wspace);
atrest.k = 100   ## restraint constant in kcal/mol/rad
atrest.Dist_ij = 8;  ## straint to a distance of 8 A, the distance between
atrest.Atom_i = 1;   ## Atom 1
atrest.Atom_j = 60;  ## and Atom 60
atrest.info()
ff.add(atrest);

## print energies as summary
ff.printEnergySummary()
ff.printEnergyByAtom();
## also show parameters
ff.info()

tra1 = OutTra_NAMD("output",wspace);
wspace.addTra(tra1)

## do some minimisation and a little MD 
min = Minimisation(ff)
min.Steps = 25
min.UpdateScr  = 1 
min.UpdateTra  = 0
min.run()

ff.printEnergySummary()
md = MolecularDynamics(ff)
md.Steps = 1000 
md.UpdateScr = 10
md.UpdateTra = 10
md.UpdateNList = 10
md.Integrator = MolecularDynamics.Langevin
md.setTargetTemp(300)
md.run()


ff.printEnergySummary()

