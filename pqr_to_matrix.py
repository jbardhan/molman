#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

parser = argparse.ArgumentParser(description = "This program takes a .pqr file (MEAD format only for now, meaning no chain field!) and writes a CRG and PDB file from it.", prog = sys.argv[0])

parser.add_argument('--pqr', metavar = 'lysozyme.pqr')
parser.add_argument('--output', metavar = '<base for output filename>')
args = parser.parse_args(sys.argv[1:])

system = mymm.Molecule()
system.read_pqr(args.pqr)

for atom in system.atoms:
    for atom2 in system.atoms:
        if atom == atom2:
            atom2.set_charge(1.0)
        else:
            atom2.set_charge(0.0)
    system.write_mead_pqr(filename=args.output+str(atom.number)+".pqr")
    
            
