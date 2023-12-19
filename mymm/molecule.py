from .atom import *
import re, os
import math
import mymm
import subprocess

class Molecule:
    def __init__(self, atoms = None, PDB = None, AMBER_CRD = None):
        self.atoms = []
        self.psf   = {}
        
        if [atoms!=None, PDB!=None, AMBER_CRD!=None].count(True) != 1:
            print("Warning: Molecule __init__ can have either atoms, PDB, or AMBER_CRD, but not more than one!\n")
            return
        
        if atoms != None:
            for atom in atoms:
                self.add_atom(atom)

        if PDB != None:
            self.read_pdb(PDB)

        if AMBER_CRD != None:
            self.read_amber_crd(AMBER_CRD)

    def get_charges(self):
        print("WARNING: reminder - Molecule.get_charges() may not return the charges in the expected order. Fix this!")
        q_vec = []
        for atom in self.atoms:
            q_vec.append(atom.q)
        return q_vec

    def write_crg(self, filename):
        fh = open(filename, 'w')
        fh.write('aaaaaarrrnnnncqqqqqqqq\n')
        for atom in self.atoms:
            fh.write(atom.crg_entry())
        fh.close()

    def assign_prmtop_names(self):
        for atom in self.atoms:
            atom.atomid = self.prmtop['atom_names'][atom.number-1]

    def assign_prmtop_radii(self, scale_factor = 1.0):
        atom_rmin2_hash = {'h1':1.3870,
                                   'h2':1.2870,
                                   'h3':1.1870,
                                   'h4':1.4090,
                                   'h5':1.3590,
                                   'ha':1.4590,
                                   'hc':1.4870,
                                   'hn':0.6000,
                                   'ho':0.0000,
                                   'hp':0.6000,
                                   'hs':0.6000,
                                   'hw':0.0000,
                                   'hx':1.1000,
                                   'o' :1.6612,
                                   'oh':1.7210,
                                   'os':1.6837,
                                   'ow':1.7683,
                                   'c' :1.9080,
                                   'c1':1.9080,
                                   'c2':1.9080,
                                   'c3':1.9080,
                                   'ca':1.9080,
                                   'cc':1.9080,
                                   'cd':1.9080,
                                   'ce':1.9080,
                                   'cf':1.9080,
                                   'cg':1.9080,
                                   'ch':1.9080,
                                   'cp':1.9080,
                                   'cq':1.9080,
                                   'cu':1.9080,
                                   'cv':1.9080,
                                   'cx':1.9080,
                                   'cy':1.9080,
                                   'cz':1.9080,
                                   'n' :1.8240,
                                   'n1':1.8240,
                                   'n2':1.8240,
                                   'n3':1.8240,
                                   'n4':1.8240,
                                   'na':1.8240,
                                   'nb':1.8240,
                                   'nc':1.8240,
                                   'nd':1.8240,
                                   'ne':1.8240,
                                   'nf':1.8240,
                                   'nh':1.8240,
                                   'no':1.8240,
                                   's' :2.0000,
                                   's2':2.0000,
                                   's4':2.0000,
                                   's6':2.0000,
                                   'sx':2.0000,
                                   'sy':2.0000,
                                   'sh':2.0000,
                                   'ss':2.0000,
                                   'p2':2.1000,
                                   'p3':2.1000,
                                   'p4':2.1000,
                                   'p5':2.1000,
                                   'pb':2.1000,
                                   'pc':2.1000,
                                   'pd':2.1000,
                                   'pe':2.1000,
                                   'pf':2.1000,
                                   'px':2.1000,
                                   'py':2.1000,
                                   'f' :1.7500,
                                   'cl':1.9480,
                                   'br':2.0200,
                                   'i' :2.1500}
        for atom in self.atoms:
            atom_radius = atom_rmin2_hash[self.prmtop['amber_atom_types'][atom.number-1]]
#            print "atom type is " + self.prmtop['amber_atom_types'][atom.number-1] + ", so radius is " + str(atom_radius)
            atom.radius = scale_factor * atom_radius
            
    def assign_prmtop_res_info(self):
        for atom in self.atoms:
            atom.resnum = len([atomnum for atomnum in self.prmtop['res_pointers'] if atomnum <= atom.number])
            atom.resid  = self.prmtop['res_labels'][atom.resnum-1]
            
    def assign_prmtop_charges(self):
        for atom in self.atoms:
            atomcharge = self.prmtop['charges'][atom.number-1]
            atom.set_charge(atomcharge)

    def read_amber_crd(self, filename):
        fh = open(filename, 'r')
        first_line = fh.readline() # junk comment ****
        second_line = fh.readline() # num atoms
        num_atoms = int(second_line.rstrip().lstrip())
        if num_atoms < 0 or num_atoms > 100000:
            print("Error!  num_atoms has illegal value " + str(num_atoms) + "!\n")
            sys.exit(1)

        atom_coords = []
#        print "read_amber_crd: num_atoms = " + str(num_atoms)
        for line in fh:
            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue
            atom_coords += line_data
                    
        fh.close()
#        print "read_amber_crd: len(atom_coords) = " + str(len(atom_coords))

        if len(atom_coords) != 3 * num_atoms:
            print("read_amber_crd: Error! len(atom_coords) != 3 * num_atoms\n")
            sys.exit(1)

        for serial in range(1, num_atoms+1):
            self.atoms.append(Atom(number = serial,
                                   atomid = "CA",
                                   resid = "ALA",
                                   segid  ="A",
                                   resnum = serial,
                                   absres = serial,
                                   x = float(atom_coords[3*(serial-1)]),
                                   y = float(atom_coords[3*(serial-1)+1]),
                                   z = float(atom_coords[3*(serial-1)+2]),
                                   q = 0))

    def read_amber_prmtop(self, PRMTOP):
        fh = open(PRMTOP, 'r')
        title     = self.read_amber_prmtop_title(fh)
        pointers  = self.read_amber_prmtop_pointers(fh)
        charges   = self.read_amber_prmtop_charges(fh)
        atom_names     = self.read_amber_prmtop_atom_names(fh)
        masses    = self.read_amber_prmtop_masses(fh)
        atom_types = self.read_amber_prmtop_atom_types(fh)
        amber_atom_types = self.read_amber_prmtop_amber_atom_types(fh)
        Acoeffs    = self.read_amber_prmtop_acoeffs(fh)
        Bcoeffs    = self.read_amber_prmtop_bcoeffs(fh)
        res_labels = self.read_amber_prmtop_res_labels(fh)
        res_pointers = self.read_amber_prmtop_res_pointers(fh)
        fh.close()

        self.prmtop = {'title':title,
                       'pointers':pointers,
                       'charges':charges,
                       'atom_names':atom_names,
                       'atom_types':atom_types,
                       'amber_atom_types':amber_atom_types,
                       'masses':masses,
                       'Acoeffs':Acoeffs,
                       'Bcoeffs':Bcoeffs,
                       'res_labels':res_labels,
                       'res_pointers':res_pointers}
        
        if False:
            print("prmtop title = " + title)
            print("prmtop pointers = " + str(pointers))
            print("prmtop charges = " + str(charges))
            print("prmtop atom_names = " + str(atom_names))
            print("prmtop masses = " + str(masses))
            print("prmtop atom_types = " + str(atom_types))
            print("prmtop Acoeffs = " + str(Acoeffs))
            print("prmtop Bcoeffs = " + str(Bcoeffs))
            print("prmtop res_labels = " + str(res_labels))
            print("prmtop res_pointers = " + str(res_pointers))
        
    
    def read_amber_prmtop_res_pointers(self, fh):
        fh.seek(0)
        in_res_pointers = False
        res_pointers = []
        while (not in_res_pointers):
            line = fh.readline()
            if (line.find('%FLAG RESIDUE_POINTER')>=0):
                in_res_pointers = True

        if in_res_pointers:
            next_line = fh.readline() # just format
            while in_res_pointers:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_res_pointers = False
                    continue
                else:
                    res_pointers += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_res_pointers: Error! No RESIDUE_POINTER section found!\n")
            sys.exit(1)
            
        res_pointers = list(map(int, res_pointers))
        return res_pointers

        
    def read_amber_prmtop_res_labels(self, fh):
        fh.seek(0)
        in_reslabels = False
        res_labels = []
        while (not in_reslabels):
            line = fh.readline()
            if (line.find('%FLAG RESIDUE_LA')>=0):
                in_reslabels = True
        if in_reslabels:
            next_line = fh.readline() # ignore, just formatting
            while in_reslabels:
                line = fh.readline()
                if (line.find('%FLAG') >=0):
                    in_reslabels = False
                    continue
                else:
                    line_data = line.rstrip().lstrip()
                    res_labels += [line_data[i:i+4] for i in range(0, len(line_data), 4)]
        else:
            print("read_amber_prmtop_res_labels: Error! No RESIDUE_LABEL section found!\n")
            sys.exit(1)
            
        return res_labels
        
    def read_amber_prmtop_acoeffs(self, fh):
        fh.seek(0)
        in_acoeffs = False
        acoeffs = []
        while (not in_acoeffs):
            line = fh.readline()
            if (line.find('%FLAG LENNARD_JONES_ACOEF')>=0):
                in_acoeffs = True

        if in_acoeffs:
            next_line = fh.readline() # just format
            while in_acoeffs:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_acoeffs = False
                    continue
                else:
                    acoeffs += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_acoeffs: Error! No LENNARD_JONES_ACOEF section found!\n")
            sys.exit(1)
            
        acoeffs = list(map(float, acoeffs))
        return acoeffs

    def read_amber_prmtop_bcoeffs(self, fh):
        fh.seek(0)
        in_bcoeffs = False
        bcoeffs = []
        while (not in_bcoeffs):
            line = fh.readline()
            if (line.find('%FLAG LENNARD_JONES_BCOEF')>=0):
                in_bcoeffs = True

        if in_bcoeffs:
            next_line = fh.readline() # just format
            while in_bcoeffs:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_bcoeffs = False
                    continue
                else:
                    bcoeffs += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_bcoeffs: Error! No LENNARD_JONES_BCOEF section found!\n")
            sys.exit(1)
            
        bcoeffs = list(map(float, bcoeffs))
        return bcoeffs

    def read_amber_prmtop_atom_types(self, fh):
        fh.seek(0)
        in_types = False
        types = []
        while (not in_types):
            line = fh.readline()
            if (line.find('%FLAG ATOM_TYPE_INDEX')>=0):
                in_types = True

        if in_types:
            next_line = fh.readline() # just format
            while in_types:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_types = False
                    continue
                else:
                    types += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_atom_types: Error! No ATOM_TYPE_INDEX section found!\n")
            sys.exit(1)
            
        types = list(map(int, types))
        return types

    def read_amber_prmtop_amber_atom_types(self, fh):
        fh.seek(0)
        in_types = False
        types = []
        while (not in_types):
            line = fh.readline()
            if (line.find('%FLAG AMBER_ATOM_TYPE')>=0):
                in_types = True

        if in_types:
            next_line = fh.readline() # just format
            while in_types:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_types = False
                    continue
                else:
                    types += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_atom_types: Error! No AMBER_ATOM_TYPE section found!\n")
            sys.exit(1)

        return types

    def read_amber_prmtop_masses(self, fh):
        fh.seek(0)
        in_masses = False
        masses = []
        while (not in_masses):
            line = fh.readline()
            if (line.find('%FLAG MASS')>=0):
                in_masses = True

        if in_masses:
            next_line = fh.readline() # just format
            while in_masses:
                line = fh.readline()
                if (line.find('%FLAG')>=0):
                    in_masses = False
                    continue
                else:
                    masses += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_masses: Error! No MASS section found!\n")
            sys.exit(1)
            
        masses = list(map(float, masses))
        return masses
    
    def read_amber_prmtop_atom_names(self, fh):
        fh.seek(0)
        in_atom_names = False
        atom_names = []
        while (not in_atom_names):
            line = fh.readline()
            if (line.find('%FLAG ATOM_NAME')>=0):
                in_atom_names = True

        if in_atom_names:
            next_line = fh.readline()  # just FORMAT
            while in_atom_names:
                line = fh.readline()
                if (line.find('%FLAG') >=0):
                    in_atom_names = False
                    continue
                else:
                    atom_names += line.rstrip().lstrip().split()
                    
        else:
            print("read_amber_prmtop_atom_names: Error! No ATOM_NAME section found!\n")
            sys.exit(1)
            
        return atom_names
        
    def read_amber_prmtop_charges(self, fh):
        charge_conversion = 18.2223
        fh.seek(0)
        in_charges = False
        charges = []
        while (not in_charges):
            line = fh.readline()
            if (line.find('%FLAG CHARGE') >= 0):
                in_charges = True

        if in_charges:
            next_line = fh.readline() # just FORMAT
            while in_charges:
                line = fh.readline()
                if (line.find('%FLAG') >= 0):
                    in_charges = False
                    continue
                else:
                    line_data = line.rstrip().lstrip().split()
                    charges += line_data
                    
        else:
            print("read_amber_prmtop_charges: Error! No CHARGE section found!\n")
            sys.exit(1)

        
        charges = list(map(float, charges))
        charges = [x/ charge_conversion for x in charges]
#        print "charges = " + str(charges)
        return charges
        
    def read_amber_prmtop_title(self, fh):
        fh.seek(0)
        in_title = False
        while (not in_title):
            line = fh.readline()
            if (line.find('%FLAG TITLE') >= 0):
                in_title = True

        if in_title:
            next_line = fh.readline()  # just %FORMAT(20a4), useless to us
            title_line = fh.readline()
            title = title_line.rstrip().lstrip()
        else:
            print("read_amber_prmtop_title: Error! No TITLE section found!\n")
            sys.exit(1)
            
        return title
            
    def read_amber_prmtop_pointers(self, fh):
        fh.seek(0)
        in_pointers = False
        pointers_keys = ['NATOM', 'NTYPES', 'NBONH', 'MBONA', 'NTHETH', 'MTHETA',
                         'NPHIH', 'MPHIA', 'NHPARM', 'NPARM', 'NNB', 'NRES',
                         'NBONA', 'NTHETA', 'NPHIA', 'NUMBND', 'NUMANG', 'NPTRA',
                         'NATYP', 'NPHB', 'IFPERT', 'NBPER', 'NGPER', 'NDPER',
                         'MBPER', 'MGPER', 'MDPER', 'IFBOX', 'NMXRS', 'IFCAP',
                         'NUMEXTRA', 'NCOPY']
        pointers_counts = []

        while (not in_pointers):
            line = fh.readline()
            if (line.find('%FLAG POINTERS') >= 0):
                in_pointers = True

        if in_pointers:
            next_line = fh.readline() # just formatting data, simple to us
            pointers_counts += fh.readline().rstrip().lstrip().split()
            pointers_counts += fh.readline().rstrip().lstrip().split()
            pointers_counts += fh.readline().rstrip().lstrip().split()
            pointers_counts += fh.readline().rstrip().lstrip().split()
        else:
            print("read_amber_prmtop_pointers: Error! No POINTERS section found!\n")
            sys.exit(1)
            
        if (len(pointers_counts) != 31) and (len(pointers_counts) != 32):
            print("Error: there should be either 31 or 32 pointers!\n")
            print("Instead there are " + str(len(pointers_counts)))
            sys.exit(1)
            
        if len(pointers_counts) == 31:
            pointers_keys.pop() # delete last element, NCOPY
            
        pointers_ints = list(map(int, pointers_counts))
        pointers_dictionary = dict(list(zip(pointers_keys, pointers_ints)))
        return pointers_dictionary
        
    def read_crd(self, filename):
        num_supposed_atoms = 0
        num_counted_atoms  = 0
        fh = open(filename, 'r')
        for line in fh:
            if (line.find('*') ==0):
                continue

            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 1:
                num_supposed_atoms = line_data[0]
                continue
            
            serial  = int(line[0:5].rstrip().lstrip())
            absres  = int(line[6:11].rstrip().lstrip())
            resName = line[11:14].rstrip().lstrip()
            name    = line[16:20].rstrip().lstrip()
            x       = float(line[20:30].rstrip().lstrip())
            y       = float(line[30:40].rstrip().lstrip())
            z       = float(line[40:50].rstrip().lstrip())
            chainID = line[50:53].rstrip().lstrip()
            resSeq  = int(line[56:60].rstrip().lstrip())
            q       = float(line[61:70].rstrip().lstrip())

            self.atoms.append(Atom(number = serial,
                                   atomid = name,
                                   resid  = resName,
                                   segid = chainID,
                                   resnum = resSeq,
                                   absres = absres,
                                   x = x,
                                   y = y,
                                   z = z,
                                   q = q))


    def read_pqr(self, filename):
        fh = open(filename, 'r')
        for line in fh:
            if not ((line.find('ATOM') == 0) or (line.find('HETATM') == 0)):
                continue
            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 10:
                (serial, name, resname, resnum, x, y, z, q, r) = tuple(line_data[1:])
                self.atoms.append(Atom(number = int(serial),
                                       atomid = name,
                                       resid  = resname,
                                       segid = "A",
                                       resnum = int(resnum),
                                       x = float(x),
                                       y = float(y),
                                       z = float(z),
                                       q = float(q),
                                       absres = "",
                                       occupancy = 0.0,
                                       tempfactor = 0.0,
                                       radius = float(r)))
            if len(line_data) == 11:
                (serial, name, resname, segname, resnum, x, y, z, q, r) = tuple(line_data[1:])
                self.atoms.append(Atom(number = int(serial),
                                       atomid = name,
                                       resid  = resname,
                                       segid = segname,
                                       resnum = int(resnum),
                                       x = float(x),
                                       y = float(y),
                                       z = float(z),
                                       q = float(q),
                                       absres = "",
                                       occupancy = 0.0,
                                       tempfactor = 0.0,
                                       radius = float(r)))

            if len(line_data) > 11:
                print("Error in read_pqr: PQR line has more than 11 fields...")
                print(("Guilty line is \n"+line))
                sys.exit(1)

        fh.close()

    def renumber(self):
        tree_rep_of_system = {}
        for atom in self.atoms:
            if atom.segid not in list(tree_rep_of_system.keys()):
                tree_rep_of_system[atom.segid] = {}

            if atom.resnum not in list(tree_rep_of_system[atom.segid].keys()):
                tree_rep_of_system[atom.segid][atom.resnum] = []

            tree_rep_of_system[atom.segid][atom.resnum].append(atom)
        
        counter = 1
        for segid in tree_rep_of_system:
            for resnum in tree_rep_of_system[segid]:
                for atom in tree_rep_of_system[segid][resnum]:
                    atom.set_number(counter)
                    counter = counter+1
        
    def add_atom(self, atom):
        if atom.__class__ != Atom:
            print("Error! Cannot add non-Atom to Molecule!\n")
            sys.exit(1)
        else:
            self.atoms.append(Atom(atom.number,atom.atomid,atom.resid,atom.segid,atom.resnum,atom.absres,atom.xyz[0],atom.xyz[1],atom.xyz[2],atom.q,atom.elem,atom.tempfactor, atom.occupancy))
            
    def translate(self, vec = None, x=None, y=None, z=None):
        for atom in self.atoms:
            atom.translate(vec, x, y, z)

    def add_molecule(self, molecule):
        if molecule.__class__ != Molecule:
            print("Error! Cannot add non-Molecule to Molecule!\n")
            sys.exit(1)
        
        for atom in molecule.atoms:
            self.add_atom(atom)
            
    def get_residues(self, selector):
        return list(set([x.resnum for x in self.atoms if eval(selector.string)]))

    def select_atoms(self, selection):
        sublist = list(filter(eval("lambda x: "+selection), self.atoms))
        return sublist
    
    def set_segid(self, selection, newsegid):
        sublist = self.select_atoms(selection)
        for atom in sublist:
            atom.segid = newsegid
            
    def set_temperature_factor(self, selection, value):
        sublist = self.select_atoms(selection)
        for atom in sublist:
            atom.tempfactor = value

    def get_box_info(self):
        minp = [999.999, 999.999, 999.999]
        maxp = [-999.999, -999.999, -999.999]
        for atom in self.atoms:
            if atom.xyz[0] < minp[0]:
                minp[0] = atom.xyz[0]
            if atom.xyz[1] < minp[1]:
                minp[1] = atom.xyz[1]
            if atom.xyz[2] < minp[2]:
                minp[2] = atom.xyz[2]
            if atom.xyz[0] > maxp[0]:
                maxp[0] = atom.xyz[0]
            if atom.xyz[1] > maxp[1]:
                maxp[1] = atom.xyz[1]
            if atom.xyz[2] > maxp[2]:
                maxp[2] = atom.xyz[2]
        boxsize = [maxp[0]-minp[0],maxp[1]-minp[1],maxp[2]-minp[2]]
        origin = [(maxp[0]+minp[0])/2.0,(maxp[1]+minp[1])/2.0,(maxp[2]+minp[2])/2.0]
        
        return (boxsize, origin)

    def parse_psf_line(self, line):
        words  = line.rstrip().lstrip().split()
        if len(words)<8:
            return (-1, 0,0,0,0,0,0,0)
        
        number = int(words[0])
        resnum = int(words[2])
        q      = float(words[6])
        mass   = float(words[7])
        return (number, words[1], resnum, words[3], words[4], words[5], q, mass)

    def read_psf(self, filename):
        fh = open(filename, 'r')
        remarks = []
        atoms = []
        modes = {'NTITLE':1, 'NATOM':2, 'NBOND':3, 'NTHETA':4, 'NPHI':5, 'NIMPHI':6,'NDON':7,'NACC':8}
        mode = 0
        for line in fh:
            if re.search('!N', line) is not None:
                for mode_name, mode_index in modes.items():
                    if re.search(mode_name, line) is not None:
                        mode = mode_index
                continue

            if mode == modes['NTITLE']:
                remarks.append(line)
            elif mode == modes['NATOM']:
                number, segid, resnum, resid, atomid, atomtype, myq, mass = self.parse_psf_line(line)
                if number > -1:
                    atoms.append(Atom(number, atomid, resid, segid, resnum, resnum, x=0.0,y=0.0,z=0.0,q=myq))
            else:
                continue # all other information is thrown away for now
        self.psf['REMARKS'] = remarks
        self.psf['ATOMS'] = atoms
        fh.close()

    def read_pdb(self, filename):
        absoluteResidues = {}
        residueCount = 0
        
        fh = open(filename,'r')
        for line in fh:
            if not ((line.find('ATOM') == 0) or (line.find('HETATM') == 0)):
                continue

            serial  = int(line[6:12].rstrip().lstrip())
            name    = line[12:16].rstrip().lstrip()
            altLoc  = line[16:17].rstrip().lstrip()
            resName = line[17:21].rstrip().lstrip()
            chainID = line[21:22].rstrip().lstrip()
            resSeq  = int(line[22:26].rstrip().lstrip())
            iCode   = line[26:27].rstrip().lstrip()
            x       = float(line[30:38].rstrip().lstrip())
            y       = float(line[38:46].rstrip().lstrip())

            if len(line) <= 54:
                z   = float(line[46:].rstrip().lstrip())
            else:
                z   = float(line[46:54].rstrip().lstrip())
                
            if len(line) < 55:
                occupancy = 1.0
            elif len(line) < 61:
                occupancy = float(line[54:].rstrip().lstrip())
            else:
                occupancy = float(line[54:60].rstrip().lstrip())

            if len(line) <= 61:
                tempfactor = 0.
            elif len(line) < 67:
                tempfactor = float(line[60:].rstrip().lstrip())
            else:
                tempfactor = float(line[60:66].rstrip().lstrip())

            if len(line) < 73:
                segid = ""
            elif len(line) <= 76:
                segid = line[72:].rstrip().lstrip()
            else:
                segid = line[72:76].rstrip().lstrip()

            if len(line) < 77:
                element = ""
            elif len(line) <= 78:
                element = line[76:].rstrip().lstrip()
            else:
                element = line[76:78].rstrip().lstrip()

            absResKey = "%s%s%s" % (chainID, resSeq, iCode)
            if absResKey not in absoluteResidues:
                residueCount+=1
                absoluteResidues[absResKey] = residueCount

            self.atoms.append(Atom(number = serial,
                                   atomid = name,
                                   resid  = resName,
                                   segid = segid,
                                   resnum = resSeq,
                                   absres = absoluteResidues[absResKey],
                                   x = x,
                                   y = y,
                                   z = z,
                                   elem = element,
                                   occupancy = occupancy,
                                   tempfactor = tempfactor))

        fh.close()

    def write_xyzr(self, filename):
        fh = open(filename,'w')
        for atom in self.atoms:
            fh.write(atom.xyzr_entry())
        fh.close()
        
    def write_crd(self, filename):
        fh = open(filename, 'w')
        fh.write("*\n")
        fh.write(" %d\n"%len(self.atoms))
        for atom in self.atoms:
            fh.write(atom.crd_entry())
        fh.close()
        
    def write_pdb(self, filename):
        fh = open(filename,'w')
        for atom in self.atoms:
            fh.write(atom.pdb_entry())
        fh.write("END\n")
        fh.close()

    def write_pdb2(self, filename):
        fh = open(filename,'w')
        for atom in self.atoms:
            fh.write(atom.pdb2_entry())
        fh.write("END\n")
        fh.close()
        
#        wholePDB = Molecule()

    def write_mead_pqr(self, filename):
        fh = open(filename, 'w')
        for atom in self.atoms:
            fh.write(atom.mead_pqr_entry())
        fh.close()

    def write_accurate_pqr(self, filename):
        fh = open(filename, 'w')
        for atom in self.atoms:
            fh.write(atom.accurate_pqr_entry())
        fh.close()

    def write_apbs_pqr(self, filename):
        fh = open(filename, 'w')
        for atom in self.atoms:
            fh.write(atom.apbs_pqr_entry())
        fh.close()

    def assign_charges(self, topology=None, patches=None):
        if topology is None:
#            print('in here!')
            if self.psf != {}:
#                print('in here too!')
                for atom in self.atoms:
                    atom.set_q_from_psf(self.psf)
                return
        else:
            for atom in self.atoms:
                atom.set_q_from_topology(topology, patches)

    def assign_radii(self, radii, patches):
        if radii is None:
            return
            
        for atom in self.atoms:
            atom.set_radius(radii, patches)

    def find_atom(self, resnum, atomid, resid):
        # warning: this just returns the first matching atom!!  really
        # you want to use a selector but I haven't implemented that
        # yet
#        print "Looking for " + resid + ":" + str(resnum) + " and atomid " + str(atomid)
        
        for atom in self.atoms:
#            print "  Current atom is " + str(atom.resid) + ":" + str(atom.resnum) + " and atomid " + str(atom.atomid )
            if atom.atomid == atomid and atom.resnum == resnum and atom.resid == resid:
                return atom
        return None

    def split_by_segment(self, output_base):
        segments = {}
    
        for atom in self.atoms:
            if atom.segid not in segments:
                segments[atom.segid] = Molecule()
                
            segments[atom.segid].add_atom(atom)

        segment_pdb_dict = {}
        for segid, segment in list(segments.items()):
            segment_pdb = output_base + segid + ".pdb"
            segment.write_pdb2(segment_pdb)
            segment_pdb_dict[segid] = {"pdb":segment_pdb,"auto":"none"}

        return segment_pdb_dict

    def calculate_dipole_moment(self):
        total_dipole = [0.0, 0.0, 0.0]
        
        for atom in self.atoms:
            atom_dipole = atom.calculate_dipole_moment()
            total_dipole[0] = total_dipole[0]+ atom_dipole[0]
            total_dipole[1] = total_dipole[1]+ atom_dipole[1]
            total_dipole[2] = total_dipole[2]+ atom_dipole[2]

        return total_dipole
    
    def rotateZ(self, angle):
        # note: angle is in degrees not radians

        for atom in self.atoms:
            atom.rotateZ(angle)

    def zero_all_charges(self):
        for atom in self.atoms:
            atom.set_charge(0.0)

    def set_group_charges(self, group_defs, residue, state):
        indices = []
        if residue['group'] not in group_defs.table.keys():
            titratable_group_def = group_defs.table['GLOBAL'][residue['group']]['atom_charge_states']
        else:
            titratable_group_def = group_defs.table[residue['group']][residue['group']]['atom_charge_states']

        selection = mymm.Selector(segids = residue['segid'], resnums=[ int(residue['resnum'])] ).string
#        print("Selection string is \""+selection+"\".\n")
        sublist = self.select_atoms(selection)
        residue=mymm.Molecule(sublist)
#        residue.write_pdb2("test.pdb")
        for atom in sublist:
            charge_state_dict = {}
            if atom.atomid in titratable_group_def.keys():
                charge_state_dict = {'zero':0.0,
                                     'neutral':titratable_group_def[atom.atomid][0],
                                     'charged':titratable_group_def[atom.atomid][1]
                                 }
                atom.set_charge(float(charge_state_dict[state]))
                indices.append(atom.number)
        return indices

    def set_titratable_group_charges(self, group_defs, residue=None, state=None):
        indices = []
        if residue is not None:
            indices.append(self.set_group_charges(group_defs, residue, state))
        else:
            for residue in group_defs.list_of_residues_to_titrate:
                indices.append(self.set_group_charges(group_defs, residue, state))

        charge_vec = self.get_charges()
        return [indices, charge_vec]
        
    def build_titratable_group_model_compound(self, group_defs, residue, state, topologyFileList, radii_list):

        if type(radii_list) is not list:
            radii_list = [radii_list]

#        print("Radii list is " + str(radii_list))
        
        build_script_vmd_file = "build.vmd"
        # store cwd
        startingdir = os.getcwd()
        
#        # make directory; die if exists or fails
#        try:
#            os.mkdir(directory)
#        except:
#            print("Error in making directory!")
#            
#        # change to directory
#        os.chdir(directory)

        bare_group_pdb = "bare_group.pdb"
        capped_group  = "capped_group"
        
        # IF residue is not global
        if residue['group'] in group_defs.table.keys():

            # create selection and make sublist of atoms
            selection = mymm.Selector(segids = residue['segid'], resnums=[ int(residue['resnum'])] ).string
            sublist = self.select_atoms(selection)
            extracted_residue =mymm.Molecule(sublist)

            # write pdb of sublist
            extracted_residue.write_pdb2(bare_group_pdb)
            
            # write buildvmd script using topologyFileList and patches, specify pdb/psf filenames
            f = open(build_script_vmd_file, "w")
            f.write("mol new " + bare_group_pdb + "\n")
            f.write("package require psfgen\n")
            if type(topologyFileList) is str:
                f.write("topology " + topologyFileList + "\n")
            else:
                for topologyFile in topologyFileList:
                    f.write("topology " + topologyFile + "\n")
            f.write("segment A {\n")
            f.write("pdb " + bare_group_pdb + "\n")
            f.write("first ace\n")
            f.write("last ct3\n")
            
            f.write("}\n")
            f.write("coordpdb " + bare_group_pdb + " A\n")
            f.write("guesscoord\n")
            f.write("writepsf " + capped_group + ".psf\n")
            f.write("writepdb " +  capped_group + ".pdb\n")
            f.write("quit\n")
            f.close()
        else:
            # ELSEIF residue is global
#            print("group is " + residue['group']+".")
#            print("keys is " + ", ".join(group_defs.table.keys()))
#            print("keys keys is " + " ".join(group_defs.table[residue['group']].keys()))
#            print("error! Termini are not implemented yet")

            # write buildvmd script using ALA and the other patch, topologyFileList and specify pdb/psf filenames
            selection = mymm.Selector(segids = residue['segid'], resnums=[ int(residue['resnum'])] ).string
            sublist = self.select_atoms(selection)
            extracted_residue =mymm.Molecule(sublist)

            # write pdb of sublist
            extracted_residue.write_pdb2(bare_group_pdb)
            
            # write buildvmd script using topologyFileList and patches, specify pdb/psf filenames
            f = open(build_script_vmd_file, "w")
            f.write("mol new " + bare_group_pdb + "\n")
            f.write("package require psfgen\n")
            if type(topologyFileList) is str:
                f.write("topology " + topologyFileList + "\n")
            else:
                for topologyFile in topologyFileList:
                    f.write("topology " + topologyFile + "\n")
            f.write("segment A {\n")
            f.write("pdb " + bare_group_pdb + "\n")
            f.write("mutate " + residue['resnum'] + " ALA \n")
            if re.match("CTER", residue['group']) is not None:
                f.write("first ace\n")
            if re.match("NTER", residue['group']) is not None:
                f.write("last ct3\n")
            
            f.write("}\n")
            f.write("coordpdb " + bare_group_pdb + " A\n")
            f.write("guesscoord\n")
            f.write("writepsf " + capped_group + ".psf\n")
            f.write("writepdb " +  capped_group + ".pdb\n")
            f.write("quit\n")
            f.close()


        # run vmd, check for warnings and errors, print both, die if errors
        build_command_list = ["vmd", "-dispdev", "text", "-e", "build.vmd"]
        output = subprocess.check_output(build_command_list, text=True)
        vmd_output_file=open("vmd_output",'w')
        vmd_output_file.write(output)
        vmd_output_file.close()

        [num_errors, warnings_and_errors] = self.process_VMD_output("vmd_output")
        if len(warnings_and_errors) > 0:
            print("There were warnings and/or errors running VMD:\n")
            print("VMD Output: " + "VMD Output:".join(warnings_and_errors))
        if num_errors > 0:
            print("Dying due to errors running VMD.\n")


        # load new molecule with PDB and PSF, assign charges,
        # set_group_charges, write out pqr and crg.
        final_sys = mymm.Molecule(PDB=capped_group+".pdb")
        final_sys.read_psf(filename=capped_group+".psf")
        final_sys.assign_charges()
        final_sys.set_titratable_group_charges(group_defs, residue, state)
        final_sys.write_crg(capped_group+".crg")

        ## assign radii
        #### write patchfile
        patch_data = mymm.Patch()
        patch_data.add_patched_site(residue['segid'],residue['resnum'],"ACE")
        patch_data.add_patched_site(residue['segid'],residue['resnum'],"CT3")

        radii_data = mymm.Radii()
        for radii_file in radii_list:
#            print("Radii file = " + radii_file)
            radii_data.read_radii_file(radii_file)
        final_sys.assign_radii(radii_data, patch_data)
        
        final_sys.write_apbs_pqr(capped_group+".pqr")
#        final_sys.write_pdb2(capped_group+".pdb2")
        # chdir back to original cwd
        os.chdir(startingdir)
        
    def process_VMD_output(self, filename):
        vmd_output_file = open(filename,'r')
        num_errors = 0
        warnings_and_errors = []
        for line in vmd_output_file:
            if re.search('error', line, re.IGNORECASE) is not None:
                warnings_and_errors.append(line)
                num_errors = num_errors+1
            if re.search('warning', line, re.IGNORECASE) is not None:
                warnings_and_errors.append(line)
                
        vmd_output_file.close()
        return [num_errors, warnings_and_errors]
    
