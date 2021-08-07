#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

parser = argparse.ArgumentParser(description="This program takes a PDB file ", prog = sys.argv[0])
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--sites', metavar = "6lyz.sites")
parser.add_argument('--top', metavar = "top_all27.inp", action = 'append')
parser.add_argument('--radii', metavar = "radii.siz", action = 'append')
parser.add_argument('--patch', metavar = 'patchfile')
args = parser.parse_args(sys.argv[1:])

args = {'pdb':'built.pdb',
        'psf':'built.psf',
        'segids':['A'],
        'isolated_subdir':'isos',
        'topologies': [global_param_dir+'/top_all22_prot_cmap.inp', proj_param_dir + '/top_titr.inp'],
        'parameters': [global_param_dir+'/par_all22_prot_cmap.inp']
        }

system = mymm.Molecule(PDB = args.pdb)
system.read_psf(filename = args.psf)

all_residues = system.get_residues(selector = mymm.Selector(segids = args['segids']))

list_of_titratable_residues = ['JR1' , 'JE1', 'JD1', 'JY1', 'JK1', 'JC1', 'JH1']
all_titration_states = {'JR1':['JR1','JR2','JR3','JR4'],
                        'JE1':['JE1','JE2'],
                        'JD1':['JD1','JD2'],
                        'JY1':['JY1','JY2'],
                        'JK1':['JK1','JK2','JK3','JK4'],
                        'JC1':['JC1','JC2'],
                        'JH1':['JH1','JH2','JH3']}

titratable_residues = {residue: system.get_resname(residue) for residue in  all_residues if system.get_resname(residue) in list_of_titratable_residues}

print "Titratable residues are: " + str(titratable_resiadues)
