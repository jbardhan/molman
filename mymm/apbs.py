from .molecule import *
import re, os, string, math, subprocess
import mymm
import copy

class ElecSection:
    def __init__(self, params):
        self.section = copy.deepcopy(params)

    def updateParams(self, new_params):
        self.section.update(new_params)
        
    def set_name(self, name):
        self.section['name'] = name

    def set_mol(self, mol):
        self.section['mol'] = mol

    def add_commands(self, commands):
        self.section['commands'].append(commands)

class Apbs:
    def __init__(self):
        self.pqrList = []
        self.elecList = []
        self.analysisList = []
        self.header_comment = ""
        self.basicSolvation = ""
        self.elecParams = {'dime': 0,
                           'cglen': 0,
                           'fglen': 0,
                           'fgcent': 0,
                           'cgcent' : 0,
                           'calc-type': 'mg-auto',
                           'mol':  0,
                           'equation':'lpbe',
                           'bcfl':'sdh',
                           'ion' : [{'charge': 1,
                                     'conc': 0.000,
                                     'radius': 2.0
                                } , {
                                    'charge': -1,
                                    'conc': 0.000,
                                    'radius': 2.0
                                 }],
                           'pdie':1.0,
                           'sdie':78.54,
                           'chgm':'spl0',
                           'srfm':'mol',
                           'srad': 0.0,
                           'swin' : 0.3,
                           'sdens' : 10.0,
                           'temp' : 298.15,
                           'commands': ['calcenergy total','calcforce no']
                       }

    def parse_APBS_output_Global_net_ELEC(self, line):
        words = line.rstrip().lstrip().split()
        energy = words[5]
        return energy

    def parse_APBS_output(self, filename):
        try:
            with open(filename,'r') as outputFile:
                lines = outputFile.read().splitlines()
        except FileNotFoundError:
            msg = "Could not open file " + filename + " for reading."
            print(msg)
            return None

        Global_lines = list(filter(lambda x: re.search("Global net ELEC energy",x), lines))
        if len(Global_lines) > 1:
            print("More than one line of "+filename+" matches \"Global net ELEC energy\":")
            print("They are:\n" + "".join(Global_lines))
            print("Parsing only the first for returning.")

        if len(Global_lines) == 0:
            print("No lines match \"Global net ELEC energy.")
            return None

        energy = self.parse_APBS_output_Global_net_ELEC(Global_lines[0])
        return energy

#    def get_elecParam(self, paramname, value_if_not_found = None):
#        if paramname in self.elecParams.keys():
#            return self.elecParams[paramname]#
#
#        return value_if_not_found

#    def set_elecParam(self, paramname, value, warning_if_already_set = None):
#        if paramname in self.elecParams.keys():
#            if self.elecParams[paramname] is not None:
#                if warning_if_already_set is not None:
#                    printf("Warning: elecParams['"+paramname+"'] is already set to "+str(self.elecParams[paramname])+"\n")
#
#        self.elecParams[paramname] = value

    def add_elecSection(self, name, molIndex, commands, additional_params_dict):
        new_elecSection = ElecSection(self.elecParams)
        new_elecSection.updateParams(additional_params_dict)
        new_elecSection.set_name(name)
        new_elecSection.set_mol(molIndex)
        new_elecSection.add_commands(commands)
        self.elecList.append(new_elecSection)

    def add_analysisSection(self, commands):
        if type(commands) is list:
            for command in commands:
                self.analysisList.append(command) 
        else:
            self.analysisList.append(commands)
    
        
    def print_APBS_input(self, filename):
        input_file = open(filename, 'w')

        self.print_APBS_header(input_file)

        self.print_APBS_read_section(input_file)
        
        for elecSection in self.elecList:
            self.print_APBS_elecSection(input_file, elecSection.section)

        for analysisSection in self.analysisList:
            self.print_APBS_analysisSection(input_file, analysisSection)

        self.print_APBS_closeout(input_file)

        input_file.close()

    def print_APBS_header(self, filehandle):
        self.print_comment_line_separator(filehandle)
        filehandle.write("# " + self.header_comment + "\n")
        self.print_comment_line_separator(filehandle)

    def print_comment_line_separator(self, filehandle):
        filehandle.write("#####"*8 + "\n")

    def print_APBS_read_section(self, filehandle):
        filehandle.write("read\n")
        for pqr in self.pqrList:
            filehandle.write("\tmol pqr " + pqr + "\n")
        filehandle.write("end\n\n")

    def print_APBS_closeout(self, filehandle):
        filehandle.write("quit\n")

    def print_APBS_analysisSection(self, filehandle, section):
        self.print_comment_line_separator(filehandle)
        filehandle.write(section + "\n")
        filehandle.write("\n")

    def print_APBS_elecSection(self, filehandle, section):
        filehandle.write("elec name " + section['name'] + "\n")
        filehandle.write("\t" + section['calc-type'] + "\n")

        if type(section['dime']) is not list:
            filehandle.write("\tdime " + (str(section['dime']) + " ")*3 + "\n")
        else:
            filehandle.write("\tdime " + " ".join(str(dim) for dim in section['dime']) + "\n")
        
        
        if type(section['cglen']) is not list:
            filehandle.write("\tcglen " + (str(section['cglen'])+" ")*3 + "\n")
        else:
            filehandle.write("\tcglen " + " ".join(str(cglen) for cglen in section['cglen']) + "\n")

        if type(section['fglen']) is not list:
            filehandle.write("\tfglen " + (str(section['fglen'])+" ")*3 + "\n")
        else:
            filehandle.write("\tfglen " + " ".join(str(fglen) for fglen in section['fglen']) + "\n")

        filehandle.write("\tcgcent mol " + str(section['cgcent']) + "\n")

        filehandle.write("\tfgcent mol " + str(section['fgcent']) + "\n")

        filehandle.write("\tmol " + str(section['mol']) + "\n")

        filehandle.write("\t" + section['equation'] + "\n")

        filehandle.write("\tbcfl " + section['bcfl'] + "\n")

        for ion in section['ion']:
            filehandle.write("\tion " + "charge " + str(ion['charge']) + " conc " +  str(ion['conc']) + " radius " + str(ion['radius']) + "\n")

        filehandle.write("\tpdie " + str(section['pdie']) + "\n")

        filehandle.write("\tsdie " + str(section['sdie']) + "\n")

        filehandle.write("\tchgm " + str(section['chgm']) + "\n")

        filehandle.write("\tsrfm " + str(section['srfm']) + "\n")

        filehandle.write("\tsrad " + str(section['srad']) + "\n")

        filehandle.write("\tswin " + str(section['swin']) + "\n")

        filehandle.write("\tsdens " + str(section['sdens']) + "\n")

        filehandle.write("\ttemp " + str(section['temp']) + "\n")

        for command in section['commands']:
            filehandle.write("\t"+command+"\n")

        filehandle.write("end\n\n")
        
