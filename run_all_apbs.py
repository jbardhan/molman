#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm

def runApbs(useSlurm = False):
    ApbsOutput = "apbs.out"
    ApbsCommand = "singularity exec --bind /run/user /tahoma/emsla60376/containers/apbs/apbs_singularity_ubuntu.sif /usr/local/bin/apbs apbs.in"
    
    if useSlurm:
        # write Slurm submission script
        submitScript = 'submit.sh'
        fh = open(submitScript, 'w')
        fh.write("#!/bin/bash\n")
        fh.write("#SBATCH --account=emsla60376\n")
        fh.write("#SBATCH --partition=normal\n")
        fh.write("#SBATCH --nodes=1\n")
        fh.write("#SBATCH --ntasks-per-node=1\n")
        fh.write("-o %j.out\n")
        fh.write("-e %j.err\n")
        fh.write("#SBATCH --time=03:00\n")
        fh.write(ApbsCommand + " > " + ApbsOutput + "\n")
        fh.close()
        # submit Slurm job with sbatch
        slurmList = ['sbatch', submitScript]
        subprocess.run(slurmList)
    else:
        ApbsListForRun = ApbsCommand.split(' ')
        with open(ApbsOutput,'w') as outfile:
            subprocess.run(ApbsListForRun, stdout=outfile)


parser = argparse.ArgumentParser()
parser.add_argument("--pdb", nargs=1, help="input PDB file")
parser.add_argument("--psf", nargs=1, help="input PSF file")
parser.add_argument("--titratables_def", nargs=1, help="file defining the titratable groups, their charges, and model pKa")
parser.add_argument("--titration_list",nargs=1, help="file listing the groups to titrate")
parser.add_argument("--topology", nargs=1, help="CHARMM style topology file--use a PARSE_prot_patch file")
parser.add_argument("--radii", nargs=1, help="radii.siz file")
parser.add_argument("--patchfile", nargs=1, help="file listing key patches")
parser.add_argument("--slurm",default=False,action='store_true')

args = parser.parse_args()

system = mymm.Molecule(PDB=args.pdb[0])
system.read_psf(filename = args.psf[0])
system.assign_charges()
systemRadii = mymm.Radii(args.radii[0])
systemPatches = mymm.Patch(args.patchfile[0])
system.assign_radii(systemRadii, systemPatches)

table = mymm.Titratable(filename = args.titratables_def[0])
table.printTitratableDefinitions()
table.readTitrationList(filename=args.titration_list[0])
table.printTitrationList()
table.validateTitrationListAgainstMolecule(system)

# Run neutral protein in the starting directory
print("Got here!")
runApbs(args.slurm)
print("Now here!")

listIndex = 0
for residue in table.listOfResiduesToTitrate:
    system.zero_all_charges()
    
    for charge_state in table.chargeStateHash.keys():
        site_plus_state = str(listIndex) + "_" + str(charge_state)
        print("*****  Running APBS calls in site plus state " + site_plus_state + "\n")
        os.chdir(site_plus_state)

        os.chdir("protein")
        runApbs(args.slurm)
        os.chdir("..")

        os.chdir("model_compound")
        runApbs(args.slurm)
        os.chdir("..")
        os.chdir("..")
            
    listIndex = listIndex+1

