#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

parser = argparse.ArgumentParser(description = "This program takes a .pqr file (MEAD format only for now, meaning no chain field!) and writes a CRG and PDB file from it.", prog = sys.argv[0])

parser.add_argument('--pqr', metavar = 'lysozyme.pqr')
parser.add_argument('--xyzr', metavar = 'lysozyme.xyzr')
args = parser.parse_args(sys.argv[1:])

pqrInfile = open(args.pqr,'r')
xyzrOutfile = open(args.xyzr,'w')
for line in pqrInfile:
    line_data = line.rstrip().lstrip().split()
    (x,y,z,q, r) = line_data[5:10]
    xyzrOutfile.write("%s %s %s %s\n"%(x,y,z,r))

