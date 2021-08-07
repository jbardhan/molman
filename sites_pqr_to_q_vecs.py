#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError, e:
        if e.errno != errno.ENOENT:
            raise

#######

def write_chargefile(filename, q, append = False):
    if append == True:
        fh = open(filename, 'a+')
    else:
        fh = open(filename, 'w')
    for curcharge in q:
        fh.write('%f '%(curcharge))
    fh.write('\n')
    fh.close()

#######
    
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

# won't do multiple chains correctly yet

parser = argparse.ArgumentParser(description="This program takes .pqr and .sites files, and generates the background charge vector and individual charged sites vectors.", prog = sys.argv[0])
parser.add_argument('--pqr', metavar = 'lysozyme.pqr')
parser.add_argument('--sites', metavar = "lysozyme.sites")
parser.add_argument('--bg', metavar = "background")
parser.add_argument('--prot', metavar = "protonated_sites")
parser.add_argument('--deprot', metavar = "deprotonated_sites")
parser.add_argument('--acidic', metavar = "acidic_sites")
args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule()
protein.read_pqr(args.pqr)

sites = mymm.Sites(args.sites)
for site in sites.site_data:
    sites.set_neutral(protein)
    
q_background = protein.get_charges()
write_chargefile(args.bg, q_background)

silentremove(args.prot)
for site in sites.site_data:
    sites.set_zero_charges(protein)
    sites.set_protonated(protein, site)
    q_one_protonated = protein.get_charges()
    write_chargefile(args.prot, q_one_protonated, append=True)

silentremove(args.deprot)
for site in sites.site_data:
    sites.set_zero_charges(protein)
    sites.set_deprotonated(protein, site)
    q_one_deprotonated = protein.get_charges()
    write_chargefile(args.deprot, q_one_deprotonated, append=True)

acidic = []
for site in sites.site_data:
    acidic.append(str(int(sites.site_types[site['type']]['acidic'])))

fh = open(args.acidic,'w')
fh.write(" ".join(acidic))
fh.close()
