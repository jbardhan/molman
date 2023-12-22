from .molecule import *
import re, os, string, math, subprocess, argparse, logging
import mymm
import copy
import numpy as np

### NOTE: the conversion constant here relies on using T=298.15 for all calculations.  It should really rest on the T in the APBS calculations

class Hybrid:
    def __init__(self, titrationDetails, proteinChargeDistributions, proteinSolvationEnergies, proteinPotentialFiles,
                 neutralChargeDistribution, neutralSolvationEnergy, neutralPotentialFiles):
        self.kbT_per_e_To_kJ_per_mol = (0.02585202*298.15/300.) * 96.485332  # first quantity from APBS docs (write.html) and second from NIST tables
        self.kJ_per_mol_To_kcal_per_mol = 1.0 / 4.184 
        self.charge_state_hash = {0: "neutral",
                                  1: "charged"
                            }
        self.titrationDetails = copy.deepcopy(titrationDetails)
        self.proteinSolvationEnergies = copy.deepcopy(proteinSolvationEnergies)
        self.proteinChargeDistributions = self.loadGroupChargeDistributions(proteinChargeDistributions)
        self.proteinPotentials = self.loadGroupPotentialFiles(proteinPotentialFiles)
        self.neutralChargeDistribution = self.loadSingleChargeDistribution(neutralChargeDistribution)
        self.neutralSolvationEnergy = copy.deepcopy(neutralSolvationEnergy)
        self.neutralPotentials = self.loadSinglePotentialFiles(neutralPotentialFiles[0], neutralPotentialFiles[1])
        self.DeltaDeltaGs = []
        self.DeltaGProteins = []
        self.DeltaGModelCompounds = []
        self.num_titratable_groups = len(self.titrationDetails.list_of_residues_to_titrate)
        self.Psi = np.zeros((self.num_titratable_groups, self.num_titratable_groups))


    def loadGroupChargeDistributions(self, chargeDistributions):
        list_index = 0
        loadedChargeDistributions = copy.deepcopy(chargeDistributions)
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            for charge_state in self.charge_state_hash.keys():
                loadedChargeDistributions[list_index][charge_state]['q']=self.loadSingleChargeDistribution(loadedChargeDistributions[list_index][charge_state]['q'])
            list_index = list_index+1

        return loadedChargeDistributions

    def loadSingleChargeDistribution(self, chargeDistribution):
        return np.array(chargeDistribution)
        
        
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


    def loadSinglePotentialFiles(self, file1, file2):
        if (re.search("ref",file1)) is not None:
            [solvPotentials, refPotentials] = self.loadAtomPotentialPair(file2, file1)
        else:
            [solvPotentials, refPotentials] = self.loadAtomPotentialPair(file1, file2)
                
        return solvPotentials-refPotentials
                
        
    def loadGroupPotentialFiles(self, potentialFiles):
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
                loadedPotentialFiles["protein"][list_index][charge_state] = self.loadSinglePotentialFiles(potentialFiles[list_index]["protein"][charge_state][0], potentialFiles[list_index]["protein"][charge_state][1])
                
                
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

            deltaQ = self.proteinChargeDistributions[list_index][1]['q']-self.proteinChargeDistributions[list_index][0]['q']
            deltaPotential = self.proteinPotentials["protein"][list_index][1]-self.proteinPotentials["protein"][list_index][0]
            explicit_Delta_G_protein = self.kbT_per_e_To_kJ_per_mol * np.inner(deltaQ, self.neutralPotentials
                                                                               + 0.5 * deltaPotential)

            self.DeltaGProteins.append(explicit_Delta_G_protein)
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
        num_neutral_charges = np.prod(self.neutralChargeDistribution.shape)
        num_neutral_potentials = np.prod(self.neutralPotentials.shape)
        print("Num neutral charges = " +  str(num_neutral_charges) + " and num neutral potentials = " + str(num_neutral_potentials))
        explicit_solvation_free_energy = 0.5 * np.inner(self.neutralPotentials, self.neutralChargeDistribution) * self.kbT_per_e_To_kJ_per_mol
        print("0.5 * phi^T q = " + str(explicit_solvation_free_energy))
        print("APBS output = " + str(self.neutralSolvationEnergy))

        list_index = 0
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]

            for charge_state in self.charge_state_hash.keys():
                explicit_free_energy = 0.5 * self.kbT_per_e_To_kJ_per_mol * np.inner(self.proteinPotentials["protein"][list_index][charge_state], self.proteinChargeDistributions[list_index][charge_state]['q'])
                print("0.5 * phi^T * q = " + str(explicit_free_energy))
                print("APBS output = " + str(self.proteinSolvationEnergies[list_index]["protein"][charge_state]))
            list_index = list_index+1


        
    def writeHybridInputFile(self, filename):
        input_file = open(filename, 'w')

        self.calculateDeltaDeltaGs()
        self.calculatePsiValues()
        list_index = 0
        #num_titratable_groups = len(self.titrationDetails.list_of_residues_to_titrate)
        input_file.write(str(self.num_titratable_groups)+"\n")

        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]               

            pKa_i_model = current_residue_info['pKa_model']
            gamma_i    = current_residue_info['gamma']

            Delta_Delta_G_i = self.DeltaDeltaGs[list_index] * self.kJ_per_mol_To_kcal_per_mol

            group_header_list = [str(x) for x in [pKa_i_model, gamma_i, Delta_Delta_G_i, list_index+1]]
            input_file.write("\t".join(group_header_list)+"\n")

            
            Psi_list = []
            for j in range(list_index+1,self.num_titratable_groups):
                Psi_i_j = self.Psi[list_index][j] * self.kJ_per_mol_To_kcal_per_mol #"("+str(list_index+1)+", "+str(j)+")"
                Psi_list.append(Psi_i_j)
                
            input_file.write("  ".join([str(x) for x in Psi_list])+"\n")
            list_index = list_index+1
                
        input_file.close()

    def calculatePsiValues(self):
        #num_titratable_groups = len(self.titrationDetails.list_of_residues_to_titrate)
        #self.Psi = np.zeros((num_titratable_groups, num_titratable_groups))
        list_index = 0

        
        for residue in self.titrationDetails.list_of_residues_to_titrate:
            if residue['group'] not in self.titrationDetails.table.keys():
                current_residue_info = self.titrationDetails.table['GLOBAL'][residue['group']]
            else:
                current_residue_info = self.titrationDetails.table[residue['group']][residue['group']]               

            gamma_i = float(current_residue_info['gamma'])
            delta_q =  self.proteinChargeDistributions[list_index][1]['q']-self.proteinChargeDistributions[list_index][0]['q']
                
            for j in range(list_index+1, self.num_titratable_groups):
                print("doing (" + str(list_index) + ", " + str(j) + ")")
                residue_j = self.titrationDetails.list_of_residues_to_titrate[j]
                if residue_j['group'] not in self.titrationDetails.table.keys():
                    current_residue_j_info = self.titrationDetails.table['GLOBAL'][residue_j['group']]
                else:
                    current_residue_j_info = self.titrationDetails.table[residue_j['group']][residue_j['group']]               

                gamma_j = float(current_residue_j_info['gamma'])
                delta_potentials = (self.proteinPotentials['protein'][j][1] - self.proteinPotentials['protein'][j][0])

                self.Psi[list_index][j] =  self.kbT_per_e_To_kJ_per_mol * np.inner(delta_q, delta_potentials) / (gamma_i * gamma_j)

                #print("delta_q = " +  str(delta_q))
                #print("delta_potentials = " + str(delta_potentials))
                #print("gamma_i = " + str(gamma_i) + " and gamma_j = " + str(gamma_j))
                #                sys.exit(1)
                
            list_index = list_index + 1
