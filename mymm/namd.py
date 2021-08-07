import subprocess, mymm, re, logging

class NamdSimulation:
    def __init__(self, default_parameter_file,
                 prefix = "",
                 additional_parameters = {},
                 files_to_source = [],
                 variables = {},
                 topologies = [],
                 commands_to_run = [],
                 min_solvation_distance = 0.0,
                 parameter_files = []):
        self.prefix = prefix
        self.build_vmd_filename = self.prefix + "build.vmd"
        self.build_vmd_output_filename = self.prefix + "build.out"
        self.build_filename = self.prefix + "built"

        self.solvate_vmd_filename = self.prefix + "solvate.vmd"
        self.solvate_vmd_output_filename = self.prefix + "solvate.out"
        self.solvate_filename  = self.prefix + "solvated"

        self.fixed_filename = None
        self.fep_filename = None

        self.equil_namd_script_filename = self.prefix+"equil.namd"
        self.equil_output_filename      = self.prefix+"equil"
        self.production_namd_script_filename = self.prefix+"prod.namd"
        self.production_output_filename = self.prefix + "prod"

        self.default_parameter_file = default_parameter_file
        self.additional_parameters = additional_parameters.copy()
        self.files_to_source = list(files_to_source)
        self.commands_to_run = list(commands_to_run)
        self.min_solvation_distance = min_solvation_distance
        self.parameter_files = list(parameter_files)
        self.variables = variables.copy()
        self.topologies = list(topologies)
        self.base_filename_list = []

        self.parameters = self.load_default_parameters(self.default_parameter_file)
        self.update_parameters(self.additional_parameters)

    def get_most_recent_base(self, extension = ""):
        return self.base_filename_list[-1] + extension

    def set_most_recent_base(self, base_filename):
        self.base_filename_list.append(base_filename)

    def write_namd_file(self, filename):
        fh = open(filename,'w')

        for key,val in self.variables.items():
            fh.write("set  %s  %s\n" % (key, str(val)))

        for parameter_file in self.parameter_files:
            fh.write("parameters  %s\n" % (parameter_file))
            
        for source_file in self.files_to_source:
            fh.write("source  %s\n" % (source_file))
            
        for key,val in self.parameters.items():
            fh.write("%s\t\t%s\n" % (key, str(val)))

        for command in self.commands_to_run:
            fh.write(command + "\n")
            
        fh.close()

    def log(self):
        logging.info('NamdSimulation')
        logging.info('\tprefix = %s',self.prefix)
        
    def load_previous(self, previous_files = None):
        if previous_files is None:
            print "Error in load_previous: you must specify outputs from previous run of NAMD!\n"
            sys.exit(1);

        for parameter in ["cellBasisVector1", "cellBasisVector2", "cellBasisVector3", "cellOrigin", "temperature"]:
            if parameter in self.parameters:
                del self.parameters[parameter]
                
        previous_data = {'extendedSystem' : previous_files['extendedSystem'] }
        
        if previous_files['bin'] is True:
            previous_data['bincoordinates'] = previous_files['bincoordinates']
            previous_data['binvelocities'] = previous_files['binvels']
        else:
            previous_data['coordinates'] = previous_files['coords']
            previous_data['velocities'] = previous_files['vels']

        self.update_parameters(previous_data)
            
    def set_commands(self, command_list):
        self.commands_to_run = list(command_list)

    def update_parameters(self, new_parameters):
        self.parameters.update(new_parameters)
        
    def load_default_parameters(self, filename):
        parameters = {}
        fh = open(filename,'r')
        for line in fh:
            words = line.rstrip().lstrip().split()
            if len(words) < 2:
                break
            parameters[words[0]] = words[1]
        fh.close()

        return parameters

    def build(self, segments, patches, guesscoord = False):
        psfgen_vmd_file = open(self.build_vmd_filename,'w')
        psfgen_vmd_file.write("package require psfgen\n")
        
        for topology_file in self.topologies:
            psfgen_vmd_file.write("topology " + topology_file + "\n")
            
        mutation_commands = {} # have to do these last
        for segment_name, segment_data in segments.items():
            psfgen_vmd_file.write("segment " + segment_name + " {\n")
            for segmentKey,segmentVal in segments[segment_name].items():
                if re.search('mutate',segmentKey) is not None:
                    mutation_commands[segmentKey] = segmentVal
                else:
                    psfgen_vmd_file.write("\t" + segmentKey + " " + segmentVal + "\n")

            for mutKey, mutVal in mutation_commands.items():
                psfgen_vmd_file.write("\t" + mutKey + " " + mutVal + "\n")
                
            psfgen_vmd_file.write("}\n")
            psfgen_vmd_file.write("coordpdb "+ segments[segment_name]["pdb"] + "\n")
                    
        for patch in patches:
            psfgen_vmd_file.write("patch " + patch + "\n")
            
        if guesscoord:
            psfgen_vmd_file.write("guesscoord\n")
                            
        psfgen_vmd_file.write("writepdb " + self.build_filename + ".pdb\n")
        psfgen_vmd_file.write("writepsf " + self.build_filename + ".psf\n")
        
        psfgen_vmd_file.write("quit\n")
        psfgen_vmd_file.close()
        
        psfgen_command_list = ["vmd", "-dispdev", "text", "-e", self.build_vmd_filename]
        
        psfgen_output = subprocess.check_output(psfgen_command_list)
        
        vmd_output = open(self.build_vmd_output_filename,'w')
        vmd_output.write(psfgen_output)
        vmd_output.close()
        
        self.set_most_recent_base(self.build_filename)

    # TO DO: need to verify that the script ran without problems!
    def get_default_solvation_options(self, base_filename):
        print "get_default_solvation_options: assuming minimum %f-A padding and cubic box.\n"%(self.min_solvation_distance)
        solute = mymm.Molecule(PDB=base_filename + ".pdb")
        boxsize, origin = solute.get_box_info()
        maxLength = max(boxsize) + 2.0 * self.min_solvation_distance
        halfL = maxLength / 2.0
        minX = origin[0] - halfL
        minY = origin[1] - halfL
        minZ = origin[2] - halfL
        maxX = origin[0] + halfL
        maxY = origin[1] + halfL
        maxZ = origin[2] + halfL
        return "-minmax {{%f %f %f} {%f %f %f}}" % (minX,minY,minZ,maxX,maxY,maxZ)

    def solvate(self, base_filename = None, 
                solvate_filename = None,
                solvation_options = None,
                solvate_vmd_filename = None,
                solvate_vmd_output_filename = None):
        if base_filename is None:
            base_filename = self.get_most_recent_base()
        if solvate_filename is None:
            solvate_filename = self.solvate_filename
        if solvation_options is None:
            solvation_options = self.get_default_solvation_options(base_filename)
        if solvate_vmd_filename is None:
            solvate_vmd_filename = self.solvate_vmd_filename
        if solvate_vmd_output_filename is None:
            solvate_vmd_output_filename = self.solvate_vmd_output_filename
            
        solvate_vmd_fh = open(solvate_vmd_filename,'w')
        solvate_vmd_fh.write("package require solvate\n")
        solvate_string = "solvate %s.psf %s.pdb %s -o %s\n" % (base_filename, base_filename, solvation_options,
                                                              solvate_filename)
        solvate_vmd_fh.write(solvate_string)
        solvate_vmd_fh.write("quit")
        solvate_vmd_fh.close()

        solvate_command_list = ["vmd", "-dispdev", "text", "-e", solvate_vmd_filename]
        solvate_output = subprocess.check_output(solvate_command_list)
        
        vmd_output = open(solvate_vmd_output_filename,'w')
        vmd_output.write(solvate_output)
        vmd_output.close()
        
        self.set_most_recent_base(solvate_filename)
        return solvate_filename

       
    def define_fixed_atoms(self, fixed_atoms_selector, column_selector = "B"):
        if column_selector != "B":
            print "Error in define_fixed_atoms(): column_selector must be B for now!\n"
            return
        
        fixed_atoms = mymm.Molecule(PDB=self.get_most_recent_base(".pdb"))
        fixed_atoms.set_temperature_factor(selection = fixed_atoms_selector, value = 1.0)
        if self.fixed_filename is None:
            self.fixed_filename = self.get_most_recent_base(".fix")
            
        fixed_atoms.write_pdb2(self.fixed_filename)
        del fixed_atoms

        fixed_atoms_data = {'fixedAtoms':"on",
                            'fixedAtomsFile':self.fixed_filename,
                            'fixedAtomsCol':column_selector}
        self.update_parameters(fixed_atoms_data)


    def define_fep_transformation(self, fep_options,
                                  appearing_atoms_selector = None,
                                  disappearing_atoms_selector = None,
                                  column_selector = "B"):
        if column_selector != "B":
            print "Error in define_fep_atoms(): column_selector must be B for now!\n"
            return

        self.fep_options = fep_options

        fep_atoms = mymm.Molecule(PDB=self.get_most_recent_base(".pdb"))
        fep_atoms.set_temperature_factor(selection = '1', value = 0.0) # see, B specific!

        if appearing_atoms_selector is not None:
            fep_atoms.set_temperature_factor(selection=appearing_atoms_selector,    value =  1.0)

        if disappearing_atoms_selector is not None:
            fep_atoms.set_temperature_factor(selection=disappearing_atoms_selector, value = -1.0)

        if self.fep_filename is None:
            self.fep_filename = self.get_most_recent_base(".fep")

        fep_atoms.write_pdb2(self.fep_filename)
        del fep_atoms

        fep_data = {'alch':"on",
                    'alchType':"FEP",
                    'alchFile':self.fep_filename,
                    'alchCol':column_selector,
                    'alchOutFile':self.fep_options.output_file,
                    'alchOutFreq':self.fep_options.output_freq,
                    'alchVdwLambdaEnd':self.fep_options.alch_vdw_lambda_end,
                    'alchElecLambdaStart':self.fep_options.alch_elec_lambda_start,
                    'alchVdwShiftCoeff':self.fep_options.alch_vdw_shift_coeff,
                    'alchDecouple':self.fep_options.alch_decouple,
                    'alchEquilSteps':self.fep_options.alch_equil_steps}
        self.update_parameters(fep_data)


    def setup_periodic_boundary_conditions(self, boxsize = None, origin = None):
        if (boxsize == None and origin == None):
            whole_system = mymm.Molecule(PDB=self.get_most_recent_base(".pdb"))
            boxsize, origin = whole_system.get_box_info()
            del whole_system
        pbc_cell = {'cellBasisVector1':"%.3f 0.0 0.0"%(boxsize[0]),
                   'cellBasisVector2':"0.0 %.3f 0.0"%(boxsize[1]),
                   'cellBasisVector3':"0.0 0.0 %.3f"%(boxsize[2]),
                   'cellOrigin': ("%.3f "*3%tuple(origin)).rstrip()}
        self.update_parameters(pbc_cell)
        
    def write_equilibration_script(self, command_list = None, namd_filename = None):

        if command_list is None:
            command_list = self.fep_options.get_equilibration_commands(self.variables['temp'])
        if namd_filename is None:
            namd_filename = self.equil_namd_script_filename
            
        self.update_parameters({'structure':self.get_most_recent_base(".psf"),
                                'coordinates':self.get_most_recent_base(".pdb"),
                                'outputname':self.equil_output_filename,
                                'restartname':self.equil_output_filename})
        self.set_commands(command_list)
        self.write_namd_file(namd_filename)
            
    def write_production_script(self, template_command_list, namd_filename = None, previous_files = None,
                                production_output_filename = None):
        if namd_filename is None:
            namd_filename = self.production_namd_script_filename
        if production_output_filename is None:
            production_output_filename = self.production_output_filename
            
        if previous_files is not None:
            self.load_previous(previous_files)
        else:
            self.load_previous({'bin':True,
                                'bincoordinates':self.equil_output_filename+".coor",
                                'binvels':self.equil_output_filename+".vel",
                                'extendedSystem':self.equil_output_filename+".xsc"})
            
        self.update_parameters({'outputname':production_output_filename,
                                'restartname':production_output_filename})

        command_list = self.fep_options.get_production_commands(template_command_list)
        self.set_commands(command_list)
        self.write_namd_file(namd_filename)


