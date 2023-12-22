import re, os, sys, subprocess, argparse, logging
import mymm

# really need to make this capable of handling multi-residue patches.

class Patch:
    def __init__(self, filename = ""):
        self.data = []
        if filename == "":
            return
        fh = open(filename, 'r')
        for line in fh:
            line_data = line.rstrip().lstrip().split()
            print("length of line data is " + str(len(line_data)))
            if len(line_data)==3:
                (segid, resnum, patch) = (line_data[0],line_data[1],line_data[2])
                self.add_patched_site(segid, resnum, patch)
        fh.close()

    def add_patched_site(self,segid, resnum, patch):
        self.data.append({'segid':segid,'resnum':resnum,'patch':patch})
    
    def get_patch(self, segid, resnum):
        for pres in self.data:
            if ( (pres['segid'] == segid)
                 and
                 (pres['resnum'] == str(resnum)) ):
                return pres['patch']
        return None

    def print_patch_data(self):
        for pres in self.data:
            print("Patch " + pres['patch'] +  " at segid " + pres['segid'] + ":" + pres['resnum'] + ".")
