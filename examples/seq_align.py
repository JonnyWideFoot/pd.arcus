from pd import *

# 1) Use a 'ForcefieldParameterSet' for naming definitions

ffps = ForcefieldParameterSet()
ffps.readlib("lib/amber03aa.ff")

s1 = BioSequence(ffps)
s1.ParseSequenceString("*s-ASP-GLU-AlA-Gly-ASP-D-GLU-ALA-GLU-A*")
flushed_print("Sequence 1: '")
s1.PrintToScreen()
flushed_print("'\n")

# 2) Use the 'MoleculeNaming' library set for naming definitions

molName = MoleculeNaming()

s2 = BioSequence(molName)
s2.ParseSequenceString("GLU-AlA-G-ASP-D-AlA-GLU")
flushed_print("Sequence 2: '")
s2.PrintToScreen()
flushed_print("'\n")

# 3) The two can be aligned...

expPair = ExpSeqPair(s1,s2)
ali = SimpleAligner(expPair)
ali.Align()
expPair.PrintToScreen()
