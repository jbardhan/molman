#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm


parser = argparse.ArgumentParser()
parser.add_argument("--pdb", nargs=1, help="input PDB file")
parser.add_argument("--psf", nargs=1, help="input PSF file")
parser.add_argument("--titratables_def", nargs=1, help="file defining the titratable groups, their charges, and model pKa")
parser.add_argument("--titration_list",nargs=1, help="file listing the groups to titrate")
parser.add_argument("--topology", nargs=1, help="CHARMM style topology file--use a PARSE_prot_patch file")
parser.add_argument("--radii", nargs=1, help="radii.siz file")
parser.add_argument("--patchfile", nargs=1, help="file listing key patches")
args = parser.parse_args()

#args = parser.parse_args('--pdb arg_capped_nter_cter_PARSE.pdb --psf arg_capped_nter_cter_PARSE.psf --titratables_def ../titratable_PARSE.def --titration_list titration_list.txt --topology /home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp --radii /home/bard415/repos/parameters/radii_fixNH3.siz --patchfile /home/bard415/repos/60376-testset1/testing-writing-subdirs/patchfile'.split())

#args = { 'pdb': 'arg_capped_nter_cter_PARSE.pdb',
#         'psf': 'arg_capped_nter_cter_PARSE.psf',
#         'def_file': '../titratable_PARSE.def',
#         'titration_list' : 'titration_list.txt',
#         'topology' : '/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp',
#         'radii' : '/home/bard415/repos/parameters/radii_fixNH3.siz',
#         'patchfile' : '/home/bard415/repos/60376-testset1/testing-writing-subdirs/patchfile'
#}

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

[proteinBoxSize,origin] = system.get_box_info()
scaleFactorCG = 4
scaleFactorFG = 2
cglenParam = [scaleFactorCG * x for x in proteinBoxSize]
fglenParam = [scaleFactorFG * x for x in proteinBoxSize]

protein = mymm.Apbs()
protein.headerComment = "Junk test APBS input"
base_dir = os.getcwd()
protein_params_hash =  {'dime': 129,
                        'cglen': cglenParam,
                        'fglen': fglenParam,
                        'fgcent': 2,
                        'cgcent': 2}
protein.addAnalysisSection(["print elecEnergy solv - ref end"])
protein.addElecSection(name="solv", molIndex=1, commands = "write atompot flat solv", additionalParamsDict=protein_params_hash)
protein_params_hash.update({'sdie':protein.elecParams['pdie']})
protein.addElecSection(name="ref", molIndex=1, commands = "write atompot flat ref", additionalParamsDict=protein_params_hash)
protein_params_hash.update({'sdie':protein.elecParams['sdie']})

system.set_titratable_group_charges(table, state="neutral")
system.write_crg("neutral_protein.crg")
system.write_apbs_pqr("neutral_protein.pqr")
neutral_protein_q_vec = system.get_charges()

protein.pqrList.append("neutral_protein.pqr")
protein.pqrList.append("neutral_protein.pqr")
protein.printApbsInput("apbs.in")
protein.pqrList.pop()
protein.pqrList.pop()

system.set_titratable_group_charges(table, state="zero")
system.write_crg("titratables_zeroed_protein.crg")
system.write_apbs_pqr("titratables_zeroed_protein.pqr")

list_index = 0
top_file = "/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp"
radii_file = "/home/bard415/repos/parameters/radii_fixNH3.siz"
charge_state_hash = {0: "neutral",
                     1: "charged"
}

numTitrationSites = len(table.listOfResiduesToTitrate)
totalNumTitrationStates = 2**numTitrationSites

for i in range(totalNumTitrationStates):
    print(f"Doing state {i} with {numTitrationSites} sites")
    # convert i to binary string and then list
    binaryString = format(i,"#0"+str(numTitrationSites+2)+"b")
    #print("binary string is " + str(binaryString))
    chargeStateList = list(binaryString)
    chargeStateList = chargeStateList[2:]
    chargeStateList.reverse()
    #print("chargeStateList = " + str([str(x) for x in chargeStateList]))
    # create directory and change to it
    exhaustiveIndexDir = "exhaustive_" + str(i)
    #sys.exit(0)
    os.mkdir(exhaustiveIndexDir)
    os.chdir(exhaustiveIndexDir)

    # loop over titratable groups and set charge state appropriately
    ### IMPORTANT: note that we don't call zero_all_charges() in this script
    j_index = 0
    for residue in table.listOfResiduesToTitrate:
        print("j_index = "+str(j_index))
        print("chargeStateList[j_index] = " + str(chargeStateList[j_index]))
        system.set_titratable_group_charges(table, residue, state=charge_state_hash[int(chargeStateList[j_index])])
        j_index=j_index+1
    
    if i==2:
        sys.exit(0)

    # assign radii and write crg, apbsPqr, apbsInput file
    system.assign_radii(systemRadii, systemPatches)
    system.write_crg(exhaustiveIndexDir+".crg")
    system.write_apbs_pqr(exhaustiveIndexDir+".pqr")

    protein.pqrList.append(exhaustiveIndexDir+".pqr")
    protein.pqrList.append(os.path.join(base_dir, "neutral_protein.pqr"))
    protein.printApbsInput("apbs.in")
    protein.pqrList.pop()
    protein.pqrList.pop()

    # change back to base directory
    os.chdir("..")

