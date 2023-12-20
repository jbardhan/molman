from .molecule import *
import re, os, string, math, subprocess, argparse, logging
import mymm
import copy

class Hybrid:
    def __init__(self, titrationDetails, proteinChargeDistributions, solvationEnergy):
        self.titrationDetails = copy.deepcopy(titrationDetails)
        self.solvationEnergies = copy.deepcopy(solvationEnergy)
        
    def writeHybridInputFile(self, filename):
        input_file = open(filename, 'w')

        list_index = 0
        num_titratable_groups = len(self.titrationDetails.list_of_residues_to_titrate)
        input_file.write(str(num_titratable_groups)+"\n")

        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]               

            pKa_i_model = current_residue_info['pKa_model']
            gamma_i    = current_residue_info['gamma']

            Delta_G_protein = 0
            Delta_G_model = self.solvationEnergies[list_index]['model_compound'][1]-self.solvationEnergies[list_index]['model_compound'][0]

            Delta_Delta_G_i = Delta_G_protein - Delta_G_model

            group_header_list = [str(x) for x in [pKa_i_model, gamma_i, Delta_Delta_G_i, list_index+1]]
            input_file.write("\t".join(group_header_list)+"\n")

            Psi_list = []
            for j in range(list_index+1,num_titratable_groups):
                Psi_i_j = "("+str(list_index+1)+", "+str(j)+")"
                Psi_list.append(Psi_i_j)
                
            input_file.write("  ".join([str(x) for x in Psi_list])+"\n")
            list_index = list_index+1
                
        input_file.close()
