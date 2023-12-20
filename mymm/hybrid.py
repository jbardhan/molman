from .molecule import *
import re, os, string, math, subprocess, argparse, logging
import mymm
import copy
import numpy as np

class Hybrid:
    def __init__(self, titrationDetails, proteinChargeDistributions, proteinSolvationEnergies, proteinPotentialFiles,
                 neutralChargeDistribution, neutralSolvationEnergy, neutralPotentialFiles):
        self.charge_state_hash = {0: "neutral",
                                  1: "charged"
                            }
        self.titrationDetails = copy.deepcopy(titrationDetails)
        self.proteinSolvationEnergies = copy.deepcopy(proteinSolvationEnergies)
        self.proteinChargeDistributions = self.loadChargeDistributions(proteinChargeDistributions)
        self.proteinPotentialFiles = self.loadPotentialFiles(proteinPotentialFiles)
        self.neutralChargeDistribution = copy.deepcopy(neutralChargeDistribution)
        self.neutralSolvationEnergy = copy.deepcopy(neutralSolvationEnergy)
        self.neutralPotentialFiles = copy.deepcopy(neutralPotentialFiles)
        self.DeltaDeltaGs = []
        self.DeltaGProteins = []
        self.DeltaGModelCompounds = []

    def loadChargeDistributions(self, chargeDistributions):
        list_index = 0
        loadedChargeDistributions = copy.deepcopy(chargeDistributions)
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            for charge_state in self.charge_state_hash.keys():
                loadedChargeDistributions[list_index][charge_state]['q']=np.array(loadedChargeDistributions[list_index][charge_state]['q'])
            list_index = list_index+1

        return loadedChargeDistributions

    def loadTextVector(self, filename):
        try:
            with open(filename, 'r') as outputfile:
                lines = outputfile.read().splitlines()
        except FileNotFoundError:
            msg = "Could not open file " + filename + " for reading."
            print(msg)
            return None

        vector = []
        for line in lines:
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]
            line_data = line.rstrip().lstrip().split()
            if len(line_data)==0:
                continue
            for value in line_data:
                vector.append(float(value))
        return vector
        
    def loadAtomPotentialPair(self, file1, file2):
        vec1 = np.array(self.loadTextVector(file1))
        vec2 = np.array(self.loadTextVector(file2))
        return vec1, vec2

    def loadPotentialFiles(self, potentialFiles):
        list_index = 0
        loadedPotentialFiles = {"protein" : {} ,
                                "model_compound" : {}
                            }
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            loadedPotentialFiles["protein"][list_index]={}
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            for charge_state in self.charge_state_hash.keys():
                if (re.search("ref",potentialFiles[list_index]["protein"][charge_state][0])) is not None:
                    [solvPotentials, refPotentials] = self.loadAtomPotentialPair(potentialFiles[list_index]["protein"][charge_state][0], potentialFiles[list_index]["protein"][charge_state][1])
                else:
                    [refPotentials, solvPotentials] = self.loadAtomPotentialPair(potentialFiles[list_index]["protein"][charge_state][0], potentialFiles[list_index]["protein"][charge_state][1])
                
                loadedPotentialFiles["protein"][list_index][charge_state] = solvPotentials-refPotentials
                
            list_index = list_index+1


        return loadedPotentialFiles

    def calculateDeltaGModelCompounds(self):
        list_index = 0
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            self.DeltaGModelCompounds.append(self.proteinSolvationEnergies[list_index]['model_compound'][1] - self.proteinSolvationEnergies[list_index]['model_compound'][0])

            list_index = list_index + 1

    
    def calculateDeltaGProteins(self):
        list_index = 0
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            self.DeltaGProteins.append(0)
            list_index = list_index + 1


        
    def calculateDeltaDeltaGs(self):
        self.calculateDeltaGProteins()
        self.calculateDeltaGModelCompounds()

        list_index = 0
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            self.DeltaDeltaGs.append(self.DeltaGProteins[list_index]-self.DeltaGModelCompounds[list_index])

            list_index = list_index + 1

    def testLoadedPotentialsAndCharges(self):
        print("Testing loaded potentials and charges\n")
        
        
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
            Delta_Delta_G_i = 0.0
            self.calculateDeltaDeltaGs()
            #            self.DeltaDeltaGs[list_index] = self.DeltaGProteins[list_index]-self.DeltaGModelCompounds[list_index]

            group_header_list = [str(x) for x in [pKa_i_model, gamma_i, Delta_Delta_G_i, list_index+1]]
            input_file.write("\t".join(group_header_list)+"\n")

            Psi_list = []
            for j in range(list_index+1,num_titratable_groups):
                Psi_i_j = "("+str(list_index+1)+", "+str(j)+")"
                Psi_list.append(Psi_i_j)
                
            input_file.write("  ".join([str(x) for x in Psi_list])+"\n")
            list_index = list_index+1
                
        input_file.close()
