from .molecule import *
import re, os, string, math, subprocess
import mymm

class Apbs:
    def __init__(self):
        self.pqrList = []
        self.basic_solvation = ""
        
    def print_APBS_input(self, filename, additional_params_dict):

        # process additional parameters -- maybe have another input dict for
        # commands??
        
        # open file
        input_file = open(filename, 'w')

        # process basic_solvation with the combined params dictionary
        # (not sure how to do this pythonically)
        
        # close file
        input_file.close()
