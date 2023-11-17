#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm

args = { 'pdb': 'arg_capped_nter_cter_PARSE.pdb',
         'psf': 'arg_capped_nter_cter_PARSE.psf',
         'def_file': '../titratable_PARSE.def',
         'titration_list' : 'titration_list.txt'}

system = mymm.Molecule(PDB=args['pdb'])
system.read_psf(filename = args['psf'])
system.assign_charges()

table = mymm.Titratable(filename = args['def_file'])
table.print_titratable_definitions()

table.read_titration_list(filename=args['titration_list'])
table.print_titration_list()
table.validate_titration_list_against_molecule(system)

system.set_titratable_group_charges(table, state="neutral")
system.write_crg("neutral_protein.crg")
system.write_apbs_pqr("neutral_protein.pqr")

system.set_titratable_group_charges(table, state="zero")
system.write_crg("titratables_zeroed_protein.crg")
system.write_apbs_pqr("titratables_zeroed_protein.pqr")

list_index = 0
for residue in table.list_of_residues_to_titrate:
    system.zero_all_charges()
    system.set_titratable_group_charges(table, residue, state="neutral") # neutral
    system.write_crg(str(list_index)+"_0.crg")
    system.write_apbs_pqr(str(list_index)+"_0.pqr")
    system.set_titratable_group_charges(table, residue, state="charged") # charged
    system.write_crg(str(list_index)+"_1.crg")
    system.write_apbs_pqr(str(list_index)+"_1.pqr")
    list_index = list_index+1

