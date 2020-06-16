from pd import *
from arcus import *

#####################################
# Begin user-customisable options
#####################################
pdbFile = "example.pdb"
startRes = 5
nRes = 6
ProduceNJoinedBranches = 100
stem = "example"
randSeed = 100
# nativeSequence = "*V-(TIKANLIFANGSTQTAEFKGTFEKATSEAYAYADTLKKDNGEWTVDVADKGYTLNIKFA)-G*"
exampleFileChainID = 'A'
sequenceOverride = "*V-(TIAAAAIFANGSTQTAEFKGTFEKATSEAYAYADTLKKDNGEWTVDVADKGYTLNIKFA)-G*"
# End options
#####################################

def createffts( wspace, useBreakableFF, summary):

    ff = Forcefield(wspace);

    if( useBreakableFF ):
        bonds = FF_BreakableBonded(wspace)
        bonds.DoBonds = False
        bonds.DoAngles = False
        ff.add( bonds )
    else:
        bonds = BondedForcefield(wspace)
        bonds.DoBonds = False
        bonds.DoAngles = False
        ff.add( bonds )

    sff = SoftVDWForcefield(wspace)
    ff.add(sff);

    if( summary ):
        ff.printEnergySummary()

    return ff;


def createff(wspace, useBreakableFF, dielec, summary ):

	ff = Forcefield(wspace)

	if( useBreakableFF ):
		bonds = FF_BreakableBonded(wspace)
		ff.add( bonds )
	else:
		bonds = BondedForcefield(wspace)
		ff.add( bonds )

	gbsa = GB_Still( wspace )
	gbsa.FastMode = True
	gbsa.DielectricSolute = dielec
	ff.add( gbsa )

	sasa = SASA_LCPO(wspace)
	sasa.GlobalASP = 0.009
	ff.add( sasa )

	if( summary ):
	   ff.printEnergySummary()

	return ff


endRes = startRes + nRes - 1

# set a standard random number seed for this test.
FastRandom.getInstance().reinitialise(randSeed)

libPath = "../../../param/"

# Load our anglsest
angset = AngleSet()
angset.loadFromFile(libPath + "angleset/big.loop.angleset")
angset.info()

# Load our forcefield parameters
ffps = FFParamSet()
ffps.readLib(libPath + "amber03aa.ff")

rotLib = RotamerLibrary( ffps )
rotLib.convertLib(libPath + "rotlib/shetty/scl-B30-occ1.0-rmsd1.0-prop20.0", RotLibConvert_Shetty() )

# Load our PDB
sys = PDB_In(ffps,pdbFile)
sys.disableRebuild() # Not critical, but will otherwise be repeated - ArcusBase will do the work.
sys.CentreOnBuild = True # Not critical, moves valid particles to the COG post-build.
bioSeq = BioSequence(ffps)
bioSeq.setTo(sequenceOverride)
sys.loadExplicit(exampleFileChainID,Polypeptide,bioSeq,0)

# Make a full workspace from the imported system
wspace = WorkSpace( sys )
wspace.info()
replacement = RotBond()
wspace.setRotatableBondList( replacement ) # Override the internal dummy rotbond array

# Write what we have interpreted from the import file	
PDB_Writer(stem + ".imported.pdb").write(wspace) # pdb

# Forcefield config
useBrokenBonded = True # Flag the use of the bonded forcefield that allows bond breaks.
ffs = createffts(wspace,useBrokenBonded,False)
ff = createff(wspace,useBrokenBonded,1.0,False)

segdef = SegmentDef(wspace,startRes,endRes,SegBreakCentre)
segdef.info()

# Loop stitcher
loopBuild = Arcus(ffs,ff,angset,rotLib)
loopBuild.UpdateNList = 5
loopBuild.OutputLevel = Verbosity.Eleven # i love that flag :-D
loopBuild.SegFanFilename = libPath + "segdist.dat"
loopBuild.ProduceNJoinedBranches = ProduceNJoinedBranches
loopBuild.buildAddExplicit( segdef ) # turn off auto-detection of what to build and explicitly define it

# override refiner to split protocol ...
storageRefiner = ArcusRefine_ToTra( stem )
loopBuild.setRefiner( storageRefiner )

print "All setup complete... Rolling the primary loop. Hold onto your hat!"

loopBuild.nameMe = stem
result = loopBuild.run()

print "Complete! :-D"
