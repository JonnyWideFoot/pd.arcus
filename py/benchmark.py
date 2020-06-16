from pd import *

def createff(psystem,ffps):
  ff = Forcefield(psystem,ffps)

  bonds = BondedForcefield()
  bonds.doimproper = 0
  ff.add( bonds ) 
  
  gbsa    = GB_Forcefield()
  gbsa.fastmode = 1
  ff.add( gbsa )

  sasa  = SASA_Forcefield()
  sasa.globalASP = 0.009
  ff.add( sasa )   
  
  return ff

# ----------------------------------------------------------------
# Benchmarking Parameters
# ----------------------------------------------------------------

itterations = 5
forcefieldFilename = "lib/amber03aa.ff"
pdbFilename = "trpcage.pdb"
chain = 'A'

# ----------------------------------------------------------------
# Application Header
# ----------------------------------------------------------------

flushed_print("")
info()
flushed_print("")

# ----------------------------------------------------------------
# Benchmarking Setup
# ----------------------------------------------------------------

# Import the required forcefield parameters
ffps = ForcefieldParameterSet()
ffps.readlib(forcefieldFilename)

# Create our system-specification on which to perform benchmarks
sim = SysSpec(ffps)
sim.add( ProteinPDB(ffps,pdbFilename,chain,0)  )

# Create a particlesystem and print info about it
psystem = ParticleSystem( sim )
psystem.info()

# Create the forcefield for this system
ff = createff(psystem,ffps)
epot_summary(psystem,ff)

resetPos = PhaseSpacePoint()
resetPos.save(psystem)

# A storage container for the clocks used in the benchmarking process
clocks = []

# ----------------------------------------------------------------
# Protocol 1: Torsional Minimisation
# ----------------------------------------------------------------

name_tm = "Torsional Minimisation"

printTitle( "Benchmarking " + name_tm )

clock_tm = StatClock(name_tm)
clock_tm.Begin()
clocks.append(clock_tm)

for i in range( itterations ):
  
    # Reload the orginal system
    resetPos.load(psystem)
    
    # Perform a minimisation
    min_tm = TorsionalMinimisation(psystem,ff)
    min_tm.InitialCapFactor = 0.1 # steric considerations
    min_tm.silent = True
    min_tm.steps = 100
    min_tm.run()
    
    clock_tm.Stamp()

    if( i == itterations -1 ):
        epot_summary(psystem,ff)

# ----------------------------------------------------------------
# Protocol 2: Cartesian Minimisation
# ----------------------------------------------------------------

name_cm = "Cartesian Minimisation"

printTitle( "Benchmarking " + name_cm )

clock_cm = StatClock(name_cm)
clock_cm.Begin()
clocks.append(clock_cm)

for i in range( itterations ):
      
    # Reload the orginal system
    resetPos.load(psystem)
    
    # Perform a minimisation
    min_cm = Minimisation(psystem,ff)
    min_cm.silent = True
    min_cm.steps = 500
    min_cm.run()
    
    clock_cm.Stamp()

    if( i == itterations -1 ):
        epot_summary(psystem,ff)

# ----------------------------------------------------------------
# Sumamrise The Benchmarking Process
# ----------------------------------------------------------------

printTitle( "Benchmarking Summary:" )

for i, clock in enumerate(clocks):
    clock.ReportSeconds()

