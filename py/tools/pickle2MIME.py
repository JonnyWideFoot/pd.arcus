#!/usr/bin/python

from pd import *
import sys

if len(sys.argv) < 3:
  print "Syntax: pickle2mime infile outfile"
  sys.exit(0)


infile = sys.argv[1]
outfile = sys.argv[2]

psp = readOldPickle(infile)
psp.writeMIME(outfile)

