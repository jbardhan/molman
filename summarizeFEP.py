#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging

import mymm

logging.info('------------------------------ START -----------------------------')

parameter_directory = os.environ['HOME']+"/Parameters/"

fep_input = mymm.ChargingFepInput(sys.argv, param_dir = parameter_directory)

fep_options = mymm.FepOptions(output_file = fep_input.fep_output_filename,
                              num_initial_equil_steps = fep_input.num_initial_equil_steps,
                              num_initial_min_steps   = fep_input.num_initial_min_steps,
                              num_total_steps_per_window = fep_input.num_total_steps_per_window,
                              num_equil_steps_per_window = fep_input.num_equil_steps_per_window)

fep = mymm.NamdSimulation(default_parameter_file =  fep_input.default_parameter_file,
                          files_to_source = fep_input.files_to_source,
                          topologies = fep_input.topologies,
                          parameter_files = fep_input.parameter_files,
                          min_solvation_distance = fep_input.min_solvation_distance,
                          variables = {'temp':fep_input.temp})

logging.info('------------------------------ END -----------------------------')
