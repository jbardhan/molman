class Mutator:
    def __init__(self, mutator_filename, full_atom_filename):
        self.residues = {}
        self.full_atom_lists = {}
        self.read_full_atom_file(full_atom_filename)
        self.read_mutator_file(mutator_filename)

    def read_full_atom_file(self, filename):
        current_residue = None
        fh = open(filename, 'r')
        for line in fh:
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]

            new_residue_start = re.search("Name", line, re.IGNORECASE)
            if new_residue_start is not None:
                line_data = line.rstrip().lstrip().split()
                current_residue = line_data[-1].upper()
                print("Now handling residue " + current_residue)
                self.full_atom_lists[current_residue] = []
                continue

            line_data = line.rstrip().lstrip().split()
            for atom in line_data:
                self.full_atom_lists[current_residue].append(atom)
        
    def read_mutator_file(self, filename):
        current_residue = None
        fh = open(filename, 'r')
        line_number = 0
        for line in fh:
            line_comment = re.search("#",line)
            if line_comment is not None:
                line = line[0:line_comment.start()]
            
            new_residue_start = re.search("Name", line, re.IGNORECASE)
            if new_residue_start is not None:
                line_data = line.rstrip().lstrip().split()
                current_residue = line_data[-1].upper()
                print("Now handling residue " + current_residue)
                self.residues[current_residue] = {}
                continue

            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue
            
            if len(line_data) != 2: 
                print("Error in mutator datafile " + filename + " at line " + str(line_number))
                print("   Two atom names per line only!  You have " + str(len(line_data)) + "!")
                return
            
            self.residues[current_residue][line_data[0]] = line_data[1]
            line_number = line_number + 1

    def print_state(self):
        print("Definitions loaded in mutator:")
        print("Protonatable residues and their full (protonated) atom lists:")
        for res_name in self.full_atom_lists:
            print(res_name + ": " + " ".join(self.full_atom_lists[res_name]))
        print("\n")
        print("Alchemical hybrid residues and the mapping between atoms:")
        for res_name in self.residues:
            print(res_name + ":")
            for key, val in self.residues[res_name].items():
                print("\t " + key + " -- " + val)

    def process_molecule(self, molecule, titr_state):
        newmolecule = Molecule()
        
        for cur_atom in molecule.atoms:
            if titr_state.atom_is_in_titrating_group(self, cur_atom):
                if self.atom_is_duplicated(titr_state, cur_atom):
#                    print "Atom is duplicated: " + str(cur_atom.atomid)
                    newmolecule.add_atom(self.get_duplicate_atom(titr_state, cur_atom))
                new_resid = titr_state.segids_titrating[cur_atom.segid][str(cur_atom.resnum)]
                cur_atom.change_resid(new_resid)

        molecule.add_molecule(newmolecule)
        molecule.renumber()

    def atom_is_duplicated(self, titr_state, atom):
#        print "residues keys is " + " ".join(self.residues.keys())
#        print "residue " + str(atom.number) + atom.segid + str(atom.resnum) + " to become " + titr_state.segids_titrating[atom.segid][str(atom.resnum)]
        if titr_state.segids_titrating[atom.segid][str(atom.resnum)] not in list(self.residues.keys()):
            return 0

        new_resid = titr_state.segids_titrating[atom.segid][str(atom.resnum)]
        if atom.atomid not in list(self.residues[new_resid].keys()):
            return 0

        return 1

    def get_disappearing_atoms(self, molecule, titr_state):
        disappearing_atoms_list = []
        for cur_atom in molecule.atoms:
            if titr_state.atom_is_in_titrating_group(self, cur_atom):
#                print "Handling atom " + str(cur_atom.number)
                if cur_atom.atomid in list(self.residues[cur_atom.resid].keys()):
                    disappearing_atoms_list.append(cur_atom.number)
        return "x.number in " + str(disappearing_atoms_list)
    
    def get_appearing_atoms(self, molecule, titr_state):
        appearing_atoms_list = []
        for cur_atom in molecule.atoms:
            if titr_state.atom_is_in_titrating_group(self, cur_atom):
                if cur_atom.atomid in list(self.residues[cur_atom.resid].values()):
                    appearing_atoms_list.append(cur_atom.number)
        return "x.number in " + str(appearing_atoms_list)

    def get_duplicate_atom(self, titr_state, atom):
        new_resid = titr_state.segids_titrating[atom.segid][str(atom.resnum)]
        return Atom(atom.number,self.residues[new_resid][atom.atomid],new_resid,atom.segid,atom.resnum,atom.absres,atom.xyz[0],atom.xyz[1],atom.xyz[2],atom.q,elem=atom.elem, tempfactor=atom.tempfactor, occupancy=atom.occupancy)

