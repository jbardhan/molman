import os, sys, subprocess, argparse, re, logging

import mymm

class ChargingFepInput:
    def __init__(self, command_line_args, param_dir):
        self.parameter_directory = param_dir
        self.default_parameter_file = self.parameter_directory + "/default.namd"
        self.files_to_source = [self.parameter_directory + "/fep.tcl"]
        self.topologies = [self.parameter_directory+"/top_all22_prot_cmap.inp",self.parameter_directory+"/sphere.rtf"]
        self.parameter_files = [self.parameter_directory + "/par_all22_prot_cmap.inp"]

        self.num_initial_equil_steps = 100000
        self.num_initial_min_steps = 1000
        self.num_total_steps_per_window = 40000
        self.num_equil_steps_per_window = 10000
        self.min_solvation_distance = 15.0
        self.temp = 300.0
        self.guesscoord = True
        self.not_water = "re.match('W',x.segid) is None"
        self.fixed_atoms = self.not_water
        self.appearing_atoms = self.not_water
        self.disappearing_atoms = "0"

        self.log_file = "createFEP.log"
        self.log_level = logging.DEBUG
        
        parser = argparse.ArgumentParser(description="This program takes a PDB file, builds, solvates, and sets up charging FEP calculations", prog = command_line_args[0])
        parser.add_argument('--pdb',    metavar = "'{'A':{'pdb':'6lyz.pdb','first':'ACE','last':'CT3'}}'", action='append')
        parser.add_argument('--fepout', metavar = "<output from FEP, default = 'fepout'>", default = 'fepout')
        parser.add_argument('--patch',  metavar = "<patchname A:1 B:2>", action='append', default = [])
        parser.add_argument('--runfep', metavar = "'runFEP start end delta'", action='append', default = [])
        parser.add_argument('--inputfile', metavar = "<input file>", default = None)
        args = parser.parse_args(command_line_args[1:])
        
        if args.inputfile is None:
            pdb_list = list(args.pdb)
            self.patches = list(args.patch) # ["DISU V:505 V:555", "DISU V:514 V:538", "DISU V:530 V:551"]
            
            self.fep_output_filename = args.fepout
            
            if len(args.runfep) == 0:
                self.runfep_commands = ['runFEP 0.0 1.0 0.05']
                
        else:
            pdb_list, self.patches, self.fep_output_filename, self.runfep_commands = self.loadInputFile(args.inputfile)

        self.segments = {} 
        for pdb in pdb_list:
            tmp = eval(pdb)
            for key,val in tmp.items():
                self.segments[key] = val
        
        logging.basicConfig(filename=self.log_file, level=self.log_level)

    def loadInputFile(self,filename):
        pdb = []
        segments = {}
        runfep_commands = []
        fep_output_filename = None
        patches = []
        
        fh = open(filename,'r')
        for line in fh:
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]

            line_data = line.rstrip().lstrip().split()
            if len(line_data) == 0:
                continue
            
            if (line_data[0] == 'pdb'):
                pdb.append(" ".join(line_data[1:]))
            elif line_data[0] == 'fepout':
                fep_output_filename = line_data[1]
            elif line_data[0] == 'patch':
                patches.append(" ".join(line_data[1:]))
            elif line_data[0] == 'runfep':
                runfep_commands.append(" ".join(line_data[1:]))
                
        fh.close()
    
        return (pdb, patches, fep_output_filename, runfep_commands)

    def log(self):
        logging.info('Parameter directory is %s', self.parameter_directory)
        logging.info('Segments = %s', str(self.segments))
        logging.info('Patches = %s', str(self.patches))
        logging.info('fep_output_filename = %s', str(self.fep_output_filename))
        logging.info('runfep_commands = %s', str(self.runfep_commands))
        
