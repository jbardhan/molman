#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging

import mymm

logging.info('------------------------------ START -----------------------------')

parameter_directory = os.environ['HOME']+"/Parameters/"

fep_input = mymm.ChargingFepInput(sys.argv, param_dir = parameter_directory)
fep_input.log()

fep_options = mymm.FepOptions(output_file = fep_input.fep_output_filename,
                              num_initial_equil_steps = fep_input.num_initial_equil_steps,
                              num_initial_min_steps   = fep_input.num_initial_min_steps,
                              num_total_steps_per_window = fep_input.num_total_steps_per_window,
                              num_equil_steps_per_window = fep_input.num_equil_steps_per_window)
#fep_options.log()

fep = mymm.NamdSimulation(default_parameter_file =  fep_input.default_parameter_file,
                          files_to_source = fep_input.files_to_source,
                          topologies = fep_input.topologies,
                          parameter_files = fep_input.parameter_files,
                          min_solvation_distance = fep_input.min_solvation_distance,
                          variables = {'temp':fep_input.temp})
fep.log()

fep.build(segments = fep_input.segments, patches = fep_input.patches, guesscoord = fep_input.guesscoord)
fep.solvate()
fep.setup_periodic_boundary_conditions()
fep.define_fixed_atoms(fep_input.fixed_atoms)
fep.define_fep_transformation(fep_options, appearing_atoms_selector = fep_input.appearing_atoms, disappearing_atoms_selector = fep_input.disappearing_atoms)
fep.write_equilibration_script()
fep.write_production_script(fep_input.runfep_commands)

# add log capabilities to NamdSimulation
# move getdefaultsolvation to logging
# selector class
# read PSF into Molecule
# checks on build and solvate
# print out which atoms are fixed in position and how many 
# print out which atoms are appearing and how many 
# print out which atoms are disappearing and how many
# print initial charge on system
# print final charge on system
# print box info

# a class to read in namd files and compare them

logging.info('------------------------------ END -----------------------------')
