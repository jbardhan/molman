import re, os, sys, subprocess, argparse, logging
import mymm

class Titratable:
    def __init__(self, filename = None):
        self.table = {}

        if filename is not None:
            self.read_titratable_definitions(filename)

    def read_titratable_definitions(self, filename):
        titration_def_file = open(filename)

        current_group = ""
        current_res   = ""  ## CAN BE GLOBAL
        
        for line in titration_def_file:
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]

            line_init_group = re.search("TITR", line)
            if line_init_group is not None:
                line_init_data = line.rstrip().lstrip().split()
                current_res   = line_init_data[1]
                current_group = line_init_data[2]
                current_gamma = line_init_data[3]
                self.table[current_res] = {}
                self.table[current_res][current_group] = {'gamma':current_gamma,
                                                          'atom_charge_states': {}}
                continue

            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue

            if len(line_data) > 3:
                print("Error, line_data should only have 3 entries!\n")

            self.table[current_res][current_group]['atom_charge_states'][line_data[0]] = [line_data[1],line_data[2]]

        titration_def_file.close()

    def print_titratable_definitions(self):
        residue_list = self.table.keys()
        print("Residues that have titratable groups:\n")
        print("\t" + ', '.join(residue_list))

        for residue in residue_list:
            print("Within residue " + residue + ":\n")

            group_list = self.table[residue].keys()
            print("\t we have " + str(len(group_list)) + " titratable group(s) defined.\n")
            net_charge_state_1 = 0.0
            net_charge_state_2 = 0.0

            for group in group_list:
                print("Within group " + group + ":\n")
                print("\t we have a gamma of " + str(self.table[residue][group]['gamma']) + " and atom table is \n")

                atom_list = self.table[residue][group]['atom_charge_states'].keys()
                for atom in atom_list:

#                    print("\t\t Atom " + atom + " = " + ', '.join(self.table[residue][group]['atom_charge_states'][atom]))
                    net_charge_state_1 = net_charge_state_1 + float(self.table[residue][group]['atom_charge_states'][atom][0])
                    net_charge_state_2 = net_charge_state_2 + float(self.table[residue][group]['atom_charge_states'][atom][1])
                print("\t Net charge on state 1 = " + str(net_charge_state_1) + "\n")
                print("\t Net charge on state 2 = " + str(net_charge_state_2) + "\n")
                net_charge_state_1 = 0.0
                net_charge_state_2 = 0.0

                print("\n")
            
