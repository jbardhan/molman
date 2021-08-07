#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

        
parser = argparse.ArgumentParser(description="This program is designed for the Mobley solvation test set.  It takes as input their AMBER .crd and .prmtop, and  then outputs a .pqr file for continuum-model calculation.", prog = sys.argv[0])
parser.add_argument('--pqr', metavar = '6lyz.pqr (output)')
parser.add_argument('--xyzr', metavar = '6lyz.xyzr (output)')
parser.add_argument('--prmtop', metavar = "6lyz.prmtop")
parser.add_argument('--crd', metavar = "6lyz.crd")
parser.add_argument('--scale', metavar = "<scale factor, default 0.92>", default = "0.92")

args = parser.parse_args(sys.argv[1:])

solute = mymm.Molecule(AMBER_CRD=args.crd)

solute.read_amber_prmtop(PRMTOP=args.prmtop)

solute.assign_prmtop_names()
solute.assign_prmtop_res_info()
solute.assign_prmtop_radii(float(args.scale))
solute.assign_prmtop_charges()

solute.write_accurate_pqr(filename = args.pqr)
solute.write_xyzr(filename = args.xyzr)

