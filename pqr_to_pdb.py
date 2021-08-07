#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

        
parser = argparse.ArgumentParser(description="This program takes a MEAD PQR file (NO SEGID FIELDS!) and converts it to PDB for psfgen.", prog = sys.argv[0])
parser.add_argument('--pqr', metavar = '6lyz.pqr')
parser.add_argument('--pdb', metavar = "6lyz.pdb")
args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule()
protein.read_pqr(args.pqr)
protein.write_pdb2(args.pdb)

