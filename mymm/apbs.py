from .molecule import *
import re, os, string, math, subprocess
import mymm
import copy

class ElecSection:
    def __init__(self, params):
        self.section = copy.deepcopy(params)

    def updateParams(self, newParams):
        self.section.update(newParams)
        
    def setName(self, name):
        self.section['name'] = name

    def setMol(self, mol):
        self.section['mol'] = mol

    def addCommands(self, commands):
        self.section['commands'].append(commands)

class ApbsMetadata:
    def __init__(self, wd=None, pqr=None, apbsIn=None, apbsOut="apbs.out", temp=None, rawFilename=None, solvFilename=None):
        self.workingDirectory = wd
        self.pqrFilename = pqr
        self.apbsInputFilename = apbsIn
        self.apbsOutputFilename = apbsOut
        self.temperature = temp
        self.rawAtomPotFilename = rawFilename
        self.solvAtomPotFilename = solvFilename

class ApbsOutputData:
    def __init__(self, apbsMetadata = None):
        self.workingDirectory = None
        self.pqrFilename = None
        self.pqr = None
        self.temperature = None
        self.rawPotentials = None
        self.solvPotentials = None
        self.reacPotentials = None

        if apbsMetadata is not None:
            self.loadFromApbsMetadata(apbsMetadata)

    def loadFromApbsMetadata(self, apbsMetadata):
        self.pqr = 0

class Apbs:
    def __init__(self, filename = None, pdie=10.0, sdie=78.54, temp=298.15):
        self.pqrList = []
        self.elecList = []
        self.analysisList = []
        self.headerComment = ""
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
                           'pdie':pdie,
                           'sdie':sdie,
                           'chgm':'spl0',
                           'srfm':'mol',
                           'srad': 0.0,
                           'swin' : 0.3,
                           'sdens' : 10.0,
                           'temp' : temp,
                           'commands': ['calcenergy total','calcforce no']
                       }
        if filename is not None:
            self.loadFromFile(filename)

    def parseApbsOutputGlobalNetELEC(self, line):
        words = line.rstrip().lstrip().split()
        energy = words[5]
        return float(energy)

    def parseApbsOutput(self, apbsOutFilename, apbsInFilename=None):
        try:
            with open(apbsOutFilename,'r') as outputFile:
                lines = outputFile.read().splitlines()
        except FileNotFoundError:
            msg = "Could not open file " + apbsOutFilename + " for reading."
            print(msg)
            return None

        inputLines = None
        if apbsInFilename is not None:
            try:
                with open(apbsInFilename,'r') as inputFile:
                    inputLines = inputFile.read().splitlines()
            except FileNotFoundError:
                msg = "Could not open file " + apbsInFilename + "for reading."
                print(msg)

        GlobalNetLines = list(filter(lambda x: re.search("Global net ELEC energy",x), lines))
        if len(GlobalNetLines) > 1:
            print("More than one line of "+apbsOutFilename+" matches \"Global net ELEC energy\":")
            print("They are:\n" + "".join(GlobalNetLines))
            print("Parsing only the first for returning.")

        if len(GlobalNetLines) == 0:
            print("No lines match \"Global net ELEC energy.")
            return None
        energy = self.parseApbsOutputGlobalNetELEC(GlobalNetLines[0])

        atompotFiles = []
        if inputLines is not None:
            atompotFlatLines = list(filter(lambda x: re.search("write atompot flat",x), inputLines))
            for line in atompotFlatLines:
                words = line.rstrip().lstrip().split()
                atompotFiles.append(os.path.join(os.getcwd(),words[3]+".txt"))

        return [energy, atompotFiles]

    def addElecSection(self, name, molIndex, commands, additionalParamsDict):
        newElecSection = ElecSection(self.elecParams)
        newElecSection.updateParams(additionalParamsDict)
        newElecSection.setName(name)
        newElecSection.setMol(molIndex)
        newElecSection.addCommands(commands)
        self.elecList.append(newElecSection)

    def addAnalysisSection(self, commands):
        if type(commands) is list:
            for command in commands:
                self.analysisList.append(command) 
        else:
            self.analysisList.append(commands)
    
    def printApbsInput(self, filename):
        inputFile = open(filename, 'w')

        self.printApbsHeader(inputFile)

        self.printApbsReadSection(inputFile)
        
        for elecSection in self.elecList:
            self.printApbsElecSection(inputFile, elecSection.section)

        for analysisSection in self.analysisList:
            self.printApbsAnalysisSection(inputFile, analysisSection)

        self.printApbsCloseout(inputFile)

        inputFile.close()

    def printApbsHeader(self, filehandle):
        self.printApbsCommentLineSeparator(filehandle)
        filehandle.write("# " + self.headerComment + "\n")
        self.printApbsCommentLineSeparator(filehandle)

    def printApbsCommentLineSeparator(self, filehandle):
        filehandle.write("#####"*8 + "\n")

    def printApbsReadSection(self, filehandle):
        filehandle.write("read\n")
        for pqr in self.pqrList:
            filehandle.write("\tmol pqr " + pqr + "\n")
        filehandle.write("end\n\n")

    def printApbsCloseout(self, filehandle):
        filehandle.write("quit\n")

    def printApbsAnalysisSection(self, filehandle, section):
        self.printApbsCommentLineSeparator(filehandle)
        filehandle.write(section + "\n")
        filehandle.write("\n")

    def printApbsElecSection(self, filehandle, section):
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
        
