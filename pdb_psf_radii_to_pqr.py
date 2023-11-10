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

patch_data = mymm.Patch(args.patch)

radii_data = mymm.Radii() 
for radii_file in args.radii:
    radii_data.read_radii_file(radii_file)

if len(args.radii) > 0:
    protein.assign_radii(radii_data, patch_data)

#if args.mead_output is True:
#    protein.write_mead_pqr(filename = args.pqr)
#else:
    protein.write_apbs_pqr(filename = args.pqr)
    
#print "Radii_data is ..."
#radii_data.print_radii_file("checkroux")
#print "Patch data is ..."
#patch_data.print_patch_data()


