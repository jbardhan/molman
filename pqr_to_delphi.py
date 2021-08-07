#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

parser = argparse.ArgumentParser(description = "This program takes a .pqr file (MEAD format only for now, meaning no chain field!) and writes a CRG and PDB file from it.", prog = sys.argv[0])

parser.add_argument('--pqr', metavar = 'lysozyme.pqr')
parser.add_argument('--output', metavar = '<base for output filename>')
args = parser.parse_args(sys.argv[1:])

system = mymm.Molecule()
system.read_pqr(args.pqr)
system.write_pdb2(args.output + ".pdb")
system.write_crg(args.output + ".crg")
