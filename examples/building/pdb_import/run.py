from pd import *

ffps = FFParamSet()

ffps.readLib("amber03aa.ff");
ffps.readLib("tip3.ff");
ffps.readLib("default.alias");
ffps.readLib("default.class");

stem = "../../pdb/"

names = [
 	"bba1.pdb",       ## Correctly complains about PYA
 	"water.pdb",      ## Correctly throws an exception - MOL is not a molecule in amber03aa.ff
 	"trpzip1.pdb",    ## DIES !
	"sh3.pdb",        ## Good
 	"acyltra.pdb",    ## DIES !
	"1fsd.pdb",       ## Good
 	"villin.pdb",     ## DIES !
	"protg.pdb",      ## Good
  "proteinAZ.pdb",  ## DIES !
	"ubifrag.pdb",    ## DIES ! 
	"trpzip.pdb",     ## Good
	"trpcage.pdb",    ## Good
	"2ci2.pdb"        ## Good
]

## Mass - try loading every PDB

for pdbname in names:
	try:
		sim = PDB_In(ffps,stem + pdbname);
		sim.load();                          ## Load everything from model 1
	except:
		continue
	
	## Write System to see what we've loaded from the PDB
	sim.info();
	sim.save( OutputFile_PDB( "output_" + pdbname ) )



