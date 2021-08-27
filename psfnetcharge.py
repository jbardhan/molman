#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

        
parser = argparse.ArgumentParser(description="Foo.", prog = sys.argv[0])
parser.add_argument('--pqr', metavar = '6lyz.pqr (output)')
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--psf', metavar = "6lyz.psf")
parser.add_argument('--radii', metavar = "radii.siz", action = 'append', default = [])
parser.add_argument('--patch', metavar = 'patchfile', default = "")
args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule(PDB=args.pdb)
protein.read_psf(filename = args.psf)
protein.assign_charges()

qvec = protein.get_charges()
netcharge = 0.0;
for q in qvec:
    netcharge += q
    
print("netcharge is " + str(netcharge))

