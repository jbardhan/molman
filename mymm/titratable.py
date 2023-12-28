import re, os, sys, subprocess, argparse, logging
import mymm

class Titratable:
    '''
    Titratable objects hold:
    - the defined titratable groups (in Titratable.table)
    - a list of what groups are defined (in Titratable.listOfDefinedGroups)
    - a list of the titratable groups in a given molecule (in Titratable.listOfResiduesToTitrate)

    Usage:
    a = Titratable("titratable_PARSE.def")

    After this call,
    a.table['GLU']['GLU'] is a hash with 
    {'gamma': -1,
     'pKa_model':4.07,
     'atom_charge_states': {
        'CD': [0.55, 0.10]
        'OE1': [-0.49, -0.55]
        'OE2': [-0.06, -0.55]
        }
    }

    The 'NTER' and 'CTER' groups can be applied to any residue and are stored within
    a.table['GLOBAL']['NTER']
    a.table['GLOBAL']['CTER']
    respectively.  (This choice made other things harder.  Better approaches are welcome and 
    worth implementing.)

    At the present time, 
    a.listOfDefinedGroups is just the RESNAMES (ASP, GLU, etc) as well as NTER and CTER.
    This is a helper because you can check if a given residue number + group involves the
    actual residue (e.g. ASP) or if it involves one of the groups listed under 'GLOBAL'.

    After you call 
    a.
    '''
    def __init__(self, filename = None):
        '''
        a = Titratable(filename = None)
        '''
        self.table = {}
        self.listOfDefinedGroups = []
        self.listOfResiduesToTitrate = []
        if filename is not None:
            self.readTitratableGroupDefinitions(filename)
        self.chargeStateHash = {0: "neutral",
                                  1: "charged"
                            }

    def readTitratableGroupDefinitions(self, filename):
        '''
        readTitratableGroupDefinitions(filename)
        '''
        TitrationDefFilehandle = open(filename)
        currentGroup = ""
        currentRes   = ""  ## CAN BE GLOBAL
        
        for line in TitrationDefFilehandle:
            lineCommentSearchResult = re.search("#", line)
            if lineCommentSearchResult is not None: 
                line = line[0:lineCommentSearchResult.start()]

            lineInitializeGroupSearchResult = re.search("TITR", line)
            if lineInitializeGroupSearchResult is not None:
                lineInitializeGroupWords = line.rstrip().lstrip().split()
                currentRes   = lineInitializeGroupWords[1]
                currentGroup = lineInitializeGroupWords[2]
                currentGamma = lineInitializeGroupWords[3]
                if len(lineInitializeGroupWords) > 4: ## include pKa of model compounds
                    currentPkaModelCompound = lineInitializeGroupWords[4]

                if currentRes not in self.table.keys():
                    self.table[currentRes] = {}

                self.table[currentRes][currentGroup] = {'gamma':currentGamma,
                                                          'pKa_model':currentPkaModelCompound,
                                                          'atom_charge_states': {}}
                self.listOfDefinedGroups.append(currentGroup)
                continue

            lineInsideGroupDefinitionWords = line.rstrip().lstrip().split()
            if len(lineInsideGroupDefinitionWords) == 0:
                continue

            if len(lineInsideGroupDefinitionWords) > 3:
                print("Error, lineInsideGroupDefinitionWords should only have 3 entries!\n")

            self.table[currentRes][currentGroup]['atom_charge_states'][lineInsideGroupDefinitionWords[0]] = [lineInsideGroupDefinitionWords[1],lineInsideGroupDefinitionWords[2]]

        TitrationDefFilehandle.close()

    def printTitratableDefinitions(self):
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
            
    def readTitrationList(self, filename):
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

            self.listOfResiduesToTitrate.append({'segid':line_data[0],
                                                'resnum':line_data[1],
                                                'group':line_data[2]
                                                })

            ## exit if this group is not recognized
            if line_data[2] not in self.listOfDefinedGroups:
                print("Error in parsing list of residues to titrate!  Group " + line_data[2] + " is not in the list of defined groups: \n")
                print("\t"+' ,'.join(self.listOfDefinedGroups) + "\n")
                sys.exit(1)

        titration_list_file.close()

    def printTitrationList(self):
        '''
        printTitrationList() simply prints the list of groups to titrate, in the form

        Titratable group CTER at A:2
        '''
        print("List of residues to titrate:\n")
        for residue in self.listOfResiduesToTitrate:
            print("Titratable group " + residue['group'] + " at " + residue['segid'] + ":" + residue['resnum']+"\n")
        print("\n")

    def validateTitrationListAgainstMolecule(self, molecule):
        '''
        validateTitrationListAgainstMolecule(protein) makes sure that the Titratable object's listOfResiduesToTitrate 
        only lists residues that actually exist in the Molecule object protein.
        '''
        for residue in self.listOfResiduesToTitrate:
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
        

        
