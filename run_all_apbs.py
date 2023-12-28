#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm

def runApbs():
    ApbsOutput = "apbs.out"
    ApbsCommand = "singularity exec --bind /run/user /tahoma/emsla60376/containers/apbs/apbs_singularity_ubuntu.sif /usr/local/bin/apbs apbs.in"
    ApbsListForRun = ApbsCommand.split(' ')
    with open(ApbsOutput,'w') as outfile:
        subprocess.run(ApbsListForRun, stdout=outfile)



args = { 'pdb': 'arg_capped_nter_cter_PARSE.pdb',
         'psf': 'arg_capped_nter_cter_PARSE.psf',
         'def_file': '../titratable_PARSE.def',
         'titration_list' : 'titration_list.txt',
         'topology' : '/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp',
         'radii' : '/home/bard415/repos/parameters/radii_fixNH3.siz',
         'patchfile' : '/home/bard415/repos/60376-testset1/testing-writing-subdirs/patchfile'
}

system = mymm.Molecule(PDB=args['pdb'])
system.read_psf(filename = args['psf'])
system.assign_charges()
systemRadii = mymm.Radii(args['radii'])
systemPatches = mymm.Patch(args['patchfile'])
system.assign_radii(systemRadii, systemPatches)

table = mymm.Titratable(filename = args['def_file'])
table.printTitratableDefinitions()
table.readTitrationList(filename=args['titration_list'])
table.printTitrationList()
table.validateTitrationListAgainstMolecule(system)

# Run neutral protein in the starting directory
print("Got here!")
runApbs()
print("Now here!")

listIndex = 0
for residue in table.listOfResiduesToTitrate:
    system.zero_all_charges()
    
    for charge_state in table.chargeStateHash.keys():
        site_plus_state = str(listIndex) + "_" + str(charge_state)
        print("*****  Running APBS calls in site plus state " + site_plus_state + "\n")
        os.chdir(site_plus_state)

        os.chdir("protein")
        runApbs()
        os.chdir("..")

        os.chdir("model_compound")
        runApbs()
        os.chdir("..")
        os.chdir("..")
            
    listIndex = listIndex+1

