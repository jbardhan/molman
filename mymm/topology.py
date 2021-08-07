import re, os, sys, subprocess, argparse, logging
import mymm

class Topology:
    def __init__(self, filename = None):
        self.header = []
        self.version  = {'major':0,'minor':0}
        self.masses   = {}
        self.declared_adjoining = []
        self.default_patches = {}
        self.autogenerate = []
        self.residues = {}
        self.patches  = {}

        if filename is not None:
            self.read_topology_file(filename)

        
    def read_topology_file(self, filename):
        fh = open(filename, 'r')
        raw_lines = fh.readlines()
        fh.close()

        self.add_header(raw_lines)
        self.add_version(raw_lines)
        self.add_masses(raw_lines)
        self.add_adjoining(raw_lines)
        self.add_default_patches(raw_lines)
        self.add_autogenerate(raw_lines)
        self.add_residues(raw_lines)
        self.add_patches(raw_lines)

        
    def add_header(self, lines):
        new_header = []
        self.header.append(new_header)

    def add_version(self, lines):
        new_version = {'major':0,'minor':0}

    def add_masses(self, lines):
        return
    
    def add_adjoining(self, lines):
        new_adjoining = []
        self.declared_adjoining.append(new_adjoining)

    def add_default_patches(self, lines):
        # doesn't do anything right now, obviously
        return

    def add_autogenerate(self, lines):
        # doesn't do anything right now, obviously
        return 

    def add_residues(self, lines):
        resid = None
        for line in lines:
            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue
            
            if re.search("RESI", line_data[0]) is not None:
                (resid, formal_charge) = (line_data[1],line_data[2])
                self.residues[resid] = {}

            if re.search("ATOM", line_data[0]) is not None:
                if resid is not None:
                    (atomid, atomtype, charge) = (line_data[1],line_data[2],line_data[3])
                    self.residues[resid][atomid] = {'type':atomtype, 'charge': float(charge)}

            if re.search("PRES", line_data[0]) is not None:
                resid = None

    def add_patches(self, lines):
        presid = None
        for line in lines:
            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue

            if re.search("RESI", line_data[0]) is not None:
                presid = None
         
            if re.search("ATOM", line_data[0]) is not None:
                if presid is not None:
                    (atomid, atomtype, charge) = (line_data[1],line_data[2],line_data[3])
                    self.patches[presid][atomid] = {'type':atomtype, 'charge': float(charge)}

            if re.search("PRES", line_data[0]) is not None:
                (presid, formal_charge) = (line_data[1],line_data[2])
                self.patches[presid] = {}

    
