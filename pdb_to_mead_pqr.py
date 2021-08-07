#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

        
parser = argparse.ArgumentParser(description="This program takes a PDB file, builds, solvates, and sets up charging FEP calculations", prog = sys.argv[0])
parser.add_argument('--pqr', metavar = '6lyz.pqr')
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--top', metavar = "top_all27.inp", action = 'append', default = [])
parser.add_argument('--radii', metavar = "radii.siz", action = 'append', default = [])
parser.add_argument('--patch', metavar = 'patchfile', default = "")
args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule(PDB=args.pdb)

patch_data = mymm.Patch(args.patch)

topology_data = mymm.Topology()
for topology_file in args.top:
    topology_data.read_topology_file(topology_file)

if len(args.top) > 0:
    protein.assign_charges(topology_data, patch_data)

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


