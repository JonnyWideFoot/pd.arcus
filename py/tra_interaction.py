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

ffps = ForcefieldParameterSet()
ffps.readlib("lib/amber03aa.ff")


flushed_print("------------------------------")
flushed_print("First Torsional Minimisation")
flushed_print("------------------------------\n")

sim = SysSpec(ffps)
sim.add( ProteinPDB(ffps,"trpcage.pdb",'A',0)  )

psystem = ParticleSystem( sim )
psystem.printPDB("read-in.pdb")
psystem.info()

tra = Traj_Out("outtra",psystem)

psystem.addtra(tra)

ff = createff(psystem,ffps)

epot_summary(psystem,ff)

min1 = TorsionalMinimisation(psystem,ff)
min1.InitialCapFactor = 0.1 # steric considerations
min1.steps = 20
min1.update_scr  = 1
min1.update_tra  = 10
min1.update_mon  = 10
min1.run()

epot_summary(psystem,ff)

flushed_print("------------------------------")
flushed_print("Second Cartesian Minimisation")
flushed_print("------------------------------\n")

# try to read the system back in from file!

sim2 = SysSpec(ffps)
loadtra( sim2, "outtra.tra", "last" )

psystem2 = ParticleSystem( sim2 )
psystem2.printPDB("read-in-2.pdb")
psystem2.info()

ff2 = createff(psystem2,ffps)

epot_summary(psystem2,ff2)

min2 = Minimisation(psystem2,ff2)
min2.steps = 40
min2.update_scr  = 10
min2.update_tra  = 10
min2.update_mon  = 10
min2.run()
   
epot_summary(psystem2,ff2)