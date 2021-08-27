import re, os, sys, subprocess, argparse, logging
import mymm

class Radii:
    def __init__(self, filename = None):
        self.data = {}
        if filename is not None:
            self.read_radii_file(filename)

    def read_radii_file(self, filename):
        fh = open(filename, 'r')
        junk = fh.readline()  

        for line in fh:
            atomid = line[0:6].rstrip().lstrip()
            resid  = line[6:9].rstrip().lstrip()
            radius = line[9:17].rstrip().lstrip()
            if len(resid) == 0:
                resid = "global"

            if resid not in self.data:
                self.data[resid] = {}
            if len(atomid) == 0:
                continue

            self.data[resid][atomid] = float(radius)
            
    def print_radii_file(self, filename = None):
        if filename is not None:
            fh = open(filename, 'w')
            sys.stdout = fh

        print("atom__res_radius_")
        for resid, radii in self.data.items():
            if resid == 'global':
                resid == "   "
            for atomid, radius in radii.items():
                print("%-6s%-4s%-8s" % (atomid, resid, radius))

        if filename is not None:
            sys.stdout = sys.__stdout__
            fh.close()

        
        
        
