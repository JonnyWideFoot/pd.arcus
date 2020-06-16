from pd import *

ffps = FFParamSet()

ffps.readLib("amber03aa.ff");

pdbname = "../../pdb/protg.pdb"

sim = PDB_In(ffps, pdbname);

###  Ok, so the sequence of ProteinG is MQYKLVINGKTLKGETTTKAVDAETAEKAFKQYANDNGVDGVWTYDDATKTFTVTE
### and we're going to map over it a mildly homologous sequence
###
###                                    MTTFKLIINGKTLKGEITIEAVDAAEAEKIFKQYANDNGIDGEWTYDDATKTFTVTE
###
### We'll let PD do the alignment and byild the model. It will reuse as much of the sidechain information to be obtained
### from the original PDB although ot all clashes can be resolved (have a look at the resultant PDB file)
###

bioseq = BioSequence();
bioseq.setTo("*M-(TTFKLIINGKTLKGEITIEAVDAAEAEKIFKQYANDNGIDGEWTYDDATKTFTVT)-E*");   # Map the particles defined in the PDB over this sequence!
sim.loadExplicit('A',Polypeptide,bioseq); # Load polypeptide from chain C using a sequece override
sim.info();                               # print some summarising information info
sim.save( OutputFile_PDB( "seq.pdb" ) )   # save a PDB file with the coordinates
