#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm
import tempfile

parser = argparse.ArgumentParser(description = "This program takes a .pqr file (MEAD format only for now, meaning no chain field) and writes a CRG and PDB file from it.", prog = sys.argv[0])

parser.add_argument('--pqr', metavar = 'ligand.pqr')
parser.add_argument('--radii', metavar = 'radii.siz')
parser.add_argument('--idiel', metavar = 'eps_prot', type = float, default = 1.0)
parser.add_argument('--odiel', metavar = 'eps_water', type = float, default = 80.0)
parser.add_argument('--saltconc', metavar = 'salt_concentration', type = float, default = 0.145)
parser.add_argument('--stern', metavar = 'Stern layer thickness', type = float, default = 2.0)
parser.add_argument('--usecavities',metavar = 'Use cavities?', default = True)
parser.add_argument('--proberad',metavar = 'Probe radius', type = float, default = 1.4)
parser.add_argument('--dieldens', metavar = 'Vertex density for dielectrics', type = float, default = 4.0)
parser.add_argument('--sterndens',metavar = 'Vertex density for Sterns', type = float, default = 4.0)
parser.add_argument('--output', metavar = '<base for output filename>')
args = parser.parse_args(sys.argv[1:])

# make temp dir, change into it

# copy PQR files here

# copy radii.siz file

# write delphi.prm file

pqr_list = []
for cur_pqr in pqr_list:
    pqr_to_xyzr(cur_pqr, temp_xyzr)
    pqr_to_delphi(cur_pqr, temp_base)
    mesh_surface()
    run_bem()
    process_bem_output()

# clean up, go home

system = mymm.Molecule()
system.read_pqr(args.pqr)
system.write_pdb2(args.output + ".pdb")
system.write_crg(args.output + ".crg")

