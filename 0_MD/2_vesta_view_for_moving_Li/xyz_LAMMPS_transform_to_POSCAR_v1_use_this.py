# KG 6/27/2021

import numpy as np
import math

Natoms=3
Ntimesteps=200001
everythis=100

fposcar = open("fposcar.txt", "w")
fvmd = open("fvmd.txt", "w")

with open("file.xyz", "r") as f:
  countl = 0
  lines = f.readlines()
  for m in range(0,Ntimesteps,1):
    countl = countl + 1 # Ignore the first unwanted line
    countl = countl + 1 # Ignore the second unwanted line
    for n in range(0,Natoms,1):
      if (m % everythis) == 0:
        templ = lines[countl].split()
        fvmd.write("Li %s %s %s\n" % (templ[1],templ[2],templ[3]))
        fposcar.write("%s %s %s Li\n" % (templ[1],templ[2],templ[3]))
      countl = countl + 1

f.close()
fvmd.close()
