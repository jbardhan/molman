#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm

parser = argparse.ArgumentParser(description="This program takes a PDB file of a structure and a mode vector, and writes out a PDB with the mode weights in the temperature factor field",
                                 prog = sys.argv[0])
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--mode', metavar = "mode1.txt")
parser.add_argument('--out', metavar = "6lyz_mode1.pdb")

args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule(PDB=args.pdb)

protein.set_temperature_factor("1", 0.0)

fh = open(args.mode,'r')
for line in fh:
    line_data = line.rstrip().lstrip().split()
    if len(line_data) != 2:
        print("Error: need two fields per line, an index and a weight!\n")
        
    atomNumber,atomWeight = line_data[0],line_data[1]
    protein.set_temperature_factor("x.number == " + atomNumber, float(atomWeight))
    
fh.close()

protein.write_pdb2(args.out)
