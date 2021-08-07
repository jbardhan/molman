#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
        
parser = argparse.ArgumentParser(description="This program takes a PDB file, builds, solvates, and sets up charging FEP calculations",
                                 prog = sys.argv[0])
parser.add_argument('--pqr', metavar = '6lyz.pqr')
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--charges', metavar = "charges.txt")
parser.add_argument('--radii', metavar = "radii.siz")
parser.add_argument('--patch', metavar = 'patchfile')
args = parser.parse_args(sys.argv[1:])

my_charges = mymm.Charges(args.charges)
my_radii   = mymm.Radii(args.radii)
my_patches = mymm.Patch(args.patch)

protein = mymm.Molecule(PDB=args.pdb)
protein.assign_charges_easy(my_charges)
