#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
proj_param_dir   = '/Users/jbardhan/repos/pka-contest/Parameters'

        
parser = argparse.ArgumentParser(description="This program takes a PDB file and list of residues to mutate, then creates the appropriate dual-topology PDB file", prog = sys.argv[0])
parser.add_argument('--atomlist', metavar = 'atoms-to-duplicate')
parser.add_argument('--reslist', metavar = 'atoms-to-duplicate')
parser.add_argument('--pdb', metavar = "6lyz.pdb")
parser.add_argument('--top', metavar = "top_all27.inp", action = 'append', default = [])
parser.add_argument('--out', metavar = "out.pdb")
args = parser.parse_args(sys.argv[1:])

protein = mymm.Molecule(PDB=args.pdb)

topology_data = mymm.Topology()
for topology_file in args.top:
    topology_data.read_topology_file(topology_file)

my_atomdata = []
fh = open(args.atomlist,'r')
for line in fh:
    line_data = line.rstrip().lstrip().split()
    my_atomdata.append({'resid':line_data[0],'atomid':line_data[1],'init':line_data[2],'final':line_data[3]})
fh.close()

my_resdata = []
fh = open(args.reslist,'r')
for line in fh:
    line_data = line.rstrip().lstrip().split()
    my_resdata.append({'index':int(line_data[0]),'orig':line_data[1],'new':line_data[2]})
fh.close()

for atom in protein.atoms:
    for mutation in my_resdata:
        if atom.resnum == mutation['index']:
            atom.resid = mutation['new']
            for atomdata in my_atomdata:
#                print "atomdata element is " + atomdata['resid'] + ":" + atomdata['atomid']
                if (atom.atomid == atomdata['atomid']):
#               if ((atom.resid == atomdata['resid']) and (atom.atomid == atomdata['atomid'])):

                    print("matched resid " + atom.resid + " and atomid " + atom.atomid)
                    atom.atomid = atomdata['init']
                    protein.add_atom(mymm.Atom(number=len(protein.atoms),
                                     atomid=atomdata['final'],
                                     resid =atom.resid,
                                     segid =atom.segid,
                                     resnum=atom.resnum,
                                     absres=atom.absres,
                                     x=atom.xyz[0],
                                     y=atom.xyz[1],
                                     z=atom.xyz[2],
                                     q=atom.q))

protein.write_pdb2(filename = args.out)
