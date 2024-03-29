#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm

args = { 'pdb': 'arg_capped_nter_cter_PARSE.pdb',
         'psf': 'arg_capped_nter_cter_PARSE.psf',
         'def_file': '../titratable_PARSE.def',
         'titration_list' : 'titration_list.txt',
         'topology' : '/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp',
         'radii' : '/home/bard415/repos/parameters/radii_fixNH3.siz'
}

system = mymm.Molecule(PDB=args['pdb'])
system.read_psf(filename = args['psf'])
system.assign_charges()

table = mymm.Titratable(filename = args['def_file'])
table.printTitratableDefinitions()

table.readTitrationList(filename=args['titration_list'])
table.printTitrationList()
table.validateTitrationListAgainstMolecule(system)

system.set_titratable_group_charges(table, state="neutral")
system.write_crg("neutral_protein.crg")
system.write_apbs_pqr("neutral_protein.pqr")

system.set_titratable_group_charges(table, state="zero")
system.write_crg("titratables_zeroed_protein.crg")
system.write_apbs_pqr("titratables_zeroed_protein.pqr")

list_index = 0
top_file = "/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp"
radii_file = "/home/bard415/repos/parameters/radii_fixNH3.siz"
charge_state_hash = {0: "neutral",
                     1: "charged"
}                     

for residue in table.listOfResiduesToTitrate:
    system.zero_all_charges()
    
    for charge_state in charge_state_hash.keys():
        site_plus_state = str(list_index) + "_" + str(charge_state)
        print("*****  Doing site plus state " + site_plus_state + "\n")
        os.mkdir(site_plus_state)
        os.chdir(site_plus_state)

        os.mkdir("protein")
        os.chdir("protein")
        system.set_titratable_group_charges(table, residue, state=charge_state_hash[charge_state]) # neutral
        system.write_crg(site_plus_state + ".crg")
        system.write_apbs_pqr(site_plus_state + ".pqr")
        os.chdir("..")

        os.mkdir("model_compound")
        os.chdir("model_compound")
        system.build_titratable_group_model_compound(table, residue, charge_state_hash[charge_state], top_file, radii_file)
        os.chdir("..")


        
        os.chdir("..") # out of site_plus_state
             
#    sys.exit(1)
    list_index = list_index+1

