import re, os, sys, subprocess, argparse, logging
import mymm

class Titratable:
    def __init__(self, filename = None):
        self.table = {}
        self.list_of_defined_groups = []
        self.list_of_residues_to_titrate = []
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
                if len(line_init_data) > 4: ## include pKa of model compounds
                    current_pka_model = line_init_data[4]

                if current_res not in self.table.keys():
                    self.table[current_res] = {}

                self.table[current_res][current_group] = {'gamma':current_gamma,
                                                          'pKa_model':current_pka_model,
                                                          'atom_charge_states': {}}
                self.list_of_defined_groups.append(current_group)
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
            
    def read_titration_list(self, filename):
        titration_list_file = open(filename)

        for line in titration_list_file:
            # strip comments
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]

            # split after whitespace stripping
            line_data = line.rstrip().lstrip().split()

            if len(line_data) == 0:
                continue
                
            if len(line_data) < 3 or len(line_data) > 3:
                print("Error in parsing list of residues to titrate!  Line\n"+line+"\n\t should have only 3 elements, a segid, a resnum, and a titratable groupname")
                sys.exit(1)

            self.list_of_residues_to_titrate.append({'segid':line_data[0],
                                                     'resnum':line_data[1],
                                                     'group':line_data[2]})

            ## exit if this group is not recognized
            if line_data[2] not in self.list_of_defined_groups:
                print("Error in parsing list of residues to titrate!  Group " + line_data[2] + " is not in the list of defined groups: \n")
                print("\t"+' ,'.join(self.list_of_defined_groups) + "\n")
                sys.exit(1)

        titration_list_file.close()

    def print_titration_list(self):
        print("List of residues to titrate:\n")
        for residue in self.list_of_residues_to_titrate:
            print("Titratable group " + residue['group'] + " at " + residue['segid'] + ":" + residue['resnum']+"\n")

        print("\n")

    def validate_titration_list_against_molecule(self, molecule):

        for residue in self.list_of_residues_to_titrate:
            thisresnums = []
            thisresnums.append(int(residue['resnum']))
            print("For titratable residue specified as " + residue['segid']+":"+residue['resnum']+":"+residue['group']+"\n")
            extracted_residue=mymm.Molecule(molecule.select_atoms(mymm.Selector(segids=residue['segid'], resnums=thisresnums).string))

            # if length = 0, report error and die
            print("\tNumber of atoms in extracted residue = "+str(len(extracted_residue.atoms))+"\n")
            if len(extracted_residue.atoms)==0:
                print("Zero atoms in Molecule match specification of titratable group!\n")
                sys.exit(1)
            
            # create list of all unique resnames
            residlist = []
            extracted_residue_atom_names = []
            for atom in extracted_residue.atoms:
                extracted_residue_atom_names.append(atom.atomid)
                if atom.resid not in residlist:
                    residlist.append(atom.resid)
                
            # SKIPPING FOR NOW if the list of resnames > 1 in length, report error and die (not sure how this would happen but)
            print("Number of unique resIDs in extracted residue = " + str(len(residlist)) + "\n")
            
            # if the titratable group IS NOT a global one (terminus), do the following
            if residue['group'] in self.table:

                # check the unique list from above matches our group 
                # if the one resname in the list of resnames does not match the residue group, report and die
                if re.match(residue['group'], residlist[0]) is None:
                    print("residue[group] does not match residlist[0]: " + residue['group'] + " vs " + residlist[0] + "\n")
                    sys.exit(1)

                # if we don't have all the atoms listed in the atom table, report error and die
                atom_table_names = self.table[residue['group']][residue['group']]['atom_charge_states'].keys()
                flag = 1
                if (all (atom_name in extracted_residue_atom_names for atom_name in atom_table_names)):
                    flag = 0
                if flag:
                    print("Not all atoms in the titratable group are present in the extracted residue!\n")
                    print("\tatom_table_names = " + ', '.join(atom_table_names)+"\n")
                    print("\textracted residue atom names = " + ', '.join(extracted_residue_atom_names) + "\n")
                    sys.exit(1)
                    
            
            # if the titratable group IS a global one (terminus), do the following
            # check if all the atoms in the titratable group atom table are also present in the extracted residue.  if so, we're fine, if not, report error and die
            if residue['group'] not in self.table:
                table_residue_name = 'GLOBAL'

                # if we don't have all the atoms listed in the atom table, report error and die
                atom_table_names = self.table[table_residue_name][residue['group']]['atom_charge_states'].keys()
                flag = 1
                if (all (atom_name in extracted_residue_atom_names for atom_name in atom_table_names)):
                    flag = 0
                if flag:
                    print("Not all atoms in the titratable group are present in the extracted residue!\n")
                    print("\tatom_table_names = " + ', '.join(atom_table_names)+"\n")
                    print("\textracted residue atom names = " + ', '.join(extracted_residue_atom_names) + "\n")
                    sys.exit(1)
                    
                ## the below doesn't work
                ##table_residue_name =  ist(self.table.keys())[list(self.table.values()).index(residue['group'])]  # this hack to get key from value is from https://stackoverflow.com/questions/8023306/get-key-by-value-in-dictionary
                
            print("... residue is okay!\n")
        print("All titratable groups verified against molecule.\n")
        

        
