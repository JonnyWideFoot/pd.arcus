from pd import *
from arcus import *

###################################
# user adjustable params
traFile = "example_5_6.tra.tra"
startRes = 5
nRes = 6
stem = "example_5_6_stage2"
randSeed = 100
# end user params
####################################

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

dielec = 1.0
endRes = startRes + nRes - 1

# set a standard random number seed for this test.
FastRandom.getInstance().reinitialise(randSeed)

libPath = "./lib/"

# Load our anglsest
angset = AngleSet()
angset.loadFromFile(libPath + "big.loop.angleset")
angset.info()

# Load our forcefield parameters
ffps = FFParamSet()
ffps.readLib(libPath + "amber03aa.ff")

rotLib = RotamerLibrary( ffps )
#rotLib.readLib( libPath + "lib/rotlib/shetty.rotamer" )
# Use the HIGH resolution lib for refinement
rotLib.convertLib(libPath + "scl-B30-occ1.0-rmsd0.5-prop20.0", RotLibConvert_Shetty() )

# Load our tra
sys = System(ffps)
traInput = Traj_In( traFile )
traInput.loadIntoSystem( sys, -1 )

# Make a full workspace from the imported system
wspace = WorkSpace( sys )
wspace.info()
replacement = RotBond() # Must stay in scope for the lifetime of this wspace...
wspace.setRotatableBondList( replacement ) # Override the internal dummy rotbond array
tra = wspace.addStdTra(stem) # tra

# Write what we have interpreted from the import file
PDB_Writer(stem + ".imported_S2.pdb").write(wspace) # pdb

# Forcefield config
useBrokenBonded = True # Flag the use of the bonded forcefield that allows bond breaks.
ffs = createffts(wspace,useBrokenBonded,False)
ff = createff(wspace,useBrokenBonded,dielec,False)

segdef = SegmentDef(wspace,startRes,endRes,SegBreakCentre)
segdef.info()

# Loop stitcher
loopBuild = ArcusReplay(ffs,ff,angset,rotLib,traInput)
loopBuild.UpdateNList = 5
loopBuild.OutputLevel = Verbosity.Eleven # I love that flag :-D
#loopBuild.enableComments( cmnt )
loopBuild.buildAddExplicit( segdef ) # turn off auto-detection of what to build and explicitly define it

# Override to split protocol ...
minimisingRefiner = ArcusRefine_CGMin()
minimisingRefiner.OutputLevel = Verbosity.Normal
minimisingRefiner.enableRotamerPack( rotLib )
minimisingRefiner.UpdateNList = 10
loopBuild.setRefiner( minimisingRefiner )

print "All setup complete... Rolling the primary loop. Hold onto your hat!"

result = loopBuild.run()

ff.calcEnergies()
ff.printEnergySummary()

# output the final structure to a PDB file
PDB_Writer(stem + ".best.pdb").write(wspace) # pdb
print "Complete! :-D\n" 
