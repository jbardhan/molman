from .molecule import *
import re, os, string, math, subprocess, argparse, logging
import mymm
import copy

class Hybrid:
    def __init__(self, titrationDetails, proteinChargeDistributions, proteinPotentials):
        self.titrationDetails = copy.deepcopy(titrationDetails)
        
    def writeHybridInputFile(self, filename):
        input_file = open(filename, 'w')

        num_titratable_groups = len(self.titratable.list_of_residues_to_titrate)
        input_file.write(str(num_titratable_groups))

        for i in range(num_titratable_groups):
###            current_residue_info = self.titratable[i]
            
###            pK_i_model = 0
###            gamma_i    = 0
###            Delta_Delta_G_i = 0

            group_header_list = [str(x) for x in [K_i_model, gamma_i, Delta_Delta_G_i, i+1]]
            input_file.write("\t".join(group_header_list))
            
            for j in range(num_titratable_groups, start = i+1):
###                Psi_i_j = 0
                
                input_file.write(str(Psi_i_j))
                
        input_file.close()
