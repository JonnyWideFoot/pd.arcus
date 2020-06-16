from pd import *
ctimeseed()      ## Set a random number seed using system time
timer()          ## Start a generic timer


## Set up a simple system containing a short peptide 
## and a simple forcefield contianing electrostatics
## and vdw forces

ffps = FFParamSet("amber03aa.ff")

sim = System(ffps)


mol =  NewProteinHelix(ffps,"*A-(APGIKKQ)-A*")
boundary = PeriodicBox(60,60,100)
sim.solvate_N(mol, boundary, 10)

sim.info()
sim.detail()
sim.printPDB("lots.pdb")

wspace = WorkSpace( sim )
rotbond = RotBond()
wspace.setRotatableBondList( rotbond )      
wspace.info()
wspace.printPDB()

tra = OutTra_NAMD("output",wspace)
wspace.addTra(tra)

ff = Forcefield(wspace)

bonds = BondedForcefield(wspace)
nb    = FF_NonBonded(wspace)

nb.Cutoff = 12.0
nb.InnerCutoff =  9.0
nb.VdwCutoff =  5.0
nb.VdwInnerCutoff =  4.0
nb.DDDielectric = True

ff.add( bonds )
ff.add( nb )

ff.printEnergySummary()


## Set up the evaluator. In this case this is simply a 
## single energy evaluation
## 

ene = Energy(ff)  ## simple evaluator (just the energy)

## Create a set of moves to be executed. Each will be called
## at every move. 

moveset = MoveSet(wspace)
moveset.add( MoleculeDisplacement( wspace, 1, 2, 0.1 ) )

mc = MonteCarlo(ene,moveset)
mc.Steps = 200               ## Set the number of steps
mc.UpdateScrAcc=True;        ## These set which events result in the
mc.UpdateScrRej=True;        ## printing of a line on the screen. (acceptances/rejections)
mc.UpdateTraAcc=True;        ## ditto for saving structures to the trajectory
mc.UpdateTraRej=True;

mc.run()


timer()

