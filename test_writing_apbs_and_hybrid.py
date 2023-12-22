#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm


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

table = mymm.Titratable(filename = args['def_file'])
table.print_titratable_definitions()

table.read_titration_list(filename=args['titration_list'])
table.print_titration_list()
table.validate_titration_list_against_molecule(system)

system.set_titratable_group_charges(table, state="neutral")
neutralProteinChargeDistribution = system.get_charges()
#system.write_crg("neutral_protein.crg")
#system.write_apbs_pqr("neutral_protein.pqr")


protein = mymm.Apbs()
[neutralProteinSolvationEnergy, neutralProteinPotentialFiles] = protein.parseApbsOutput("apbs.out","apbs.in")
#print("neutral_protein solvation free energy = " + str(neutralProteinSolvationEnergy) + " kJ/mol")
#print("neutral protein potential files are ")
#for newfile in neutralProteinPotentialFiles:
#        print("\t" + newfile)
#sys.exit(1)

      
protein.headerComment = "Junk test APBS input"
base_dir = os.getcwd()
protein_params_hash =  {'dime': 33,
                        'cglen': [40, 40, 40],
                        'fglen': [20, 20, 20],
                        'fgcent': 2,
                        'cgcent':2}
protein.addAnalysisSection(["print elecEnergy solv - ref end"])
protein.addElecSection(name="solv", molIndex=1, commands = "write atompot flat solv", additional_params_dict=protein_params_hash)
protein_params_hash.update({'sdie':protein.elecParams['pdie']})
protein.addElecSection(name="ref", molIndex=1, commands = "write atompot flat ref", additional_params_dict=protein_params_hash)
protein_params_hash.update({'sdie':protein.elecParams['sdie']})


system.set_titratable_group_charges(table, state="zero")
systemRadii = mymm.Radii(args['radii'])
systemPatches = mymm.Patch(args['patchfile'])

system.assign_radii(systemRadii, systemPatches)
#system.write_crg("titratables_zeroed_protein.crg")
#system.write_apbs_pqr("titratables_zeroed_protein.pqr")

list_index = 0
top_file = "/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp"
radii_file = "/home/bard415/repos/parameters/radii_fixNH3.siz"
charge_state_hash = {0: "neutral",
                     1: "charged"
}

proteinChargeDistributions = {}
solvationEnergy = {}
potentialFiles = {}
for residue in table.listOfResiduesToTitrate:
    system.zero_all_charges()

    proteinChargeDistributions[list_index] = {}
    solvationEnergy[list_index] = {'protein':{},
                                   'model_compound':{}
                               }

    potentialFiles[list_index] = {'protein':{},
                                   'model_compound':{}
                               }

    for charge_state in charge_state_hash.keys():
        site_plus_state = str(list_index) + "_" + str(charge_state)
        print("*****  Doing site plus state " + site_plus_state + "\n")
#        os.mkdir(site_plus_state)
        os.chdir(site_plus_state)

#        os.mkdir("protein")
        os.chdir("protein")
        [atom_indices, charge_vec] = system.set_titratable_group_charges(table, residue, state=charge_state_hash[charge_state]) 
        system.assign_radii(systemRadii, systemPatches)
        
 #       print("assigned radii\n")
#        system.write_crg(site_plus_state + ".crg")
#        system.write_apbs_pqr(site_plus_state + ".pqr")

        protein.pqrList.append(site_plus_state+".pqr")
        protein.pqrList.append(os.path.join(base_dir, "neutral_protein.pqr"))
#        protein.printApbsInput("apbs.in")
        protein.pqrList.pop()
        protein.pqrList.pop()

        proteinChargeDistributions[list_index][charge_state] = {'indices': atom_indices,
                                                                'q': charge_vec}
        [solvationEnergy[list_index]['protein'][charge_state],potentialFiles[list_index]['protein'][charge_state]]=protein.parseApbsOutput("apbs.out","apbs.in")

        os.chdir("..")
                     
#        os.mkdir("model_compound")
        os.chdir("model_compound") ### CONFUSING BECAUSE SYSTEM BUILDTGMC CHANGES DIRS, REMOVE THAT FUNCTIONALITY OR CREATE EQUIV ABSTRCTION FOR PROTEIN 
#        system.build_titratable_group_model_compound(table, residue, charge_state_hash[charge_state], top_file, radii_file)
        protein.pqrList.append("capped_group.pqr")
        protein.pqrList.append(os.path.join(base_dir, "neutral_protein.pqr"))
#        protein.printApbsInput("apbs.in")
        [solvationEnergy[list_index]['model_compound'][charge_state],potentialFiles[list_index]['model_compound'][charge_state]]=protein.parseApbsOutput("apbs.out","apbs.in")
        protein.pqrList.pop()
        protein.pqrList.pop()
        os.chdir("..") # out of model_compound

        os.chdir("..") # out of site_plus_state
             
#    sys.exit(1)
    list_index = list_index+1

for siteIndex in proteinChargeDistributions.keys():
    for chargeStateIndex in proteinChargeDistributions[siteIndex].keys():
        print("for " + str(siteIndex) + "_" + str(chargeStateIndex) + ":")
        indices = [str(x) for x in proteinChargeDistributions[siteIndex][chargeStateIndex]['indices']]
        print("\t the " + str(len(indices)) + " indices of the relevant group are " + " ".join(indices) )
        print("\t and the q vector for this charge state is \n" +"\n".join([str(x) for x in proteinChargeDistributions[siteIndex][chargeStateIndex]['q']]))

myhybrid = mymm.Hybrid(table, proteinChargeDistributions, solvationEnergy, potentialFiles, neutralProteinChargeDistribution, neutralProteinSolvationEnergy, neutralProteinPotentialFiles)
myhybrid.writeHybridInputFile("hybrid.out")
myhybrid.testLoadedPotentialsAndCharges()
##############################################
