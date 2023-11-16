#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
#proj_param_dir   = os.getenv('HOME') + "/repos/parameters"

args = {'pdb':'built.pdb',
        'psf':'built.psf',
        'pqr':'built.pqr',
        'group_defs':'titratable_groups.def',
        'titratable_list':'titratable_list.txt',
        'patchinfo':'patchfile',
        'segids':['X'],
        'isolated_subdir':'isos',
        'topologies': [global_param_dir+'/top_all22_prot_cmap.inp', proj_param_dir + '/top_titr.inp'],
        'parameters': [global_param_dir+'/par_all22_prot_cmap.inp']
        }

try:
    os.mkdir(args['isolated_subdir'])
except:
    pass

system = mymm.Molecule(PDB = args['pdb'])
system.read_psf(filename = args['psf'])

####
all_residues = system.get_residues(selector = mymm.Selector(segids = args['segids']))
titratable_groups_definitions = system.read_titration_group_defs(args['group_defs'])

protein_titratable_list = system.read_groups_to_titrate(args['titratable_list'])

#list_of_titratable_residues = ['JR1' , 'JE1', 'JD1', 'JY1', 'JK1', 'JC1', 'JH1']
#all_titration_states = {'JR1':['JR1','JR2','JR3','JR4'],
                        'JE1':['JE1','JE2'],
                        'JD1':['JD1','JD2'],
                        'JY1':['JY1','JY2'],
                        'JK1':['JK1','JK2','JK3','JK4'],
                        'JC1':['JC1','JC2'],
                        'JH1':['JH1','JH2','JH3']}
#titratable_residues = {residue: system.get_resname(residue) for residue in  all_residues if system.get_resname(residue) in list_of_titratable_residues}

print "Titratable residues are: " + str(titratable_resiadues)
sys.exit(0)

#####
                       
for cur_residue in titratable_residues:
    isolated_residue = mymm.Molecule(system.select_atoms(mymm.Selector(segids = args['segids'], resnums = [cur_residue]).string))

    # save current dir. make new directory. change to it
    cwd  = os.getcwd()
    dir_for_isolated_residue = args['isolated_subdir'] + '/resnum_%s'%(cur_residue)
    try:
        os.mkdir(dir_for_isolated_residue)
    except:
        pass
    os.chdir(dir_for_isolated_residue)
    
    # write PDB file, making note of residue segid
    pdb_name = 'resnum_%s_prebuild.pdb'%(cur_residue)
    isolated_residue.write_pdb2(pdb_name)
    isolated_segid = isolated_residue.atoms[0].segid

    # build using NamdSimulation
    new_isolated = mymm.NamdSimulation(default_parameter_file = global_param_dir+"/default.namd",
                                       topologies = args['topologies'],
                                       parameter_files = args['parameters'])
    new_isolated.build(segments = {isolated_segid:{'pdb':pdb_name, 'first':'ACE','last':'CT3'}},
                       patches = {},
                       guesscoord = True)
                                
    # BEM solver specific stuff
    write_and_call_setup_and_mesh(pdb = 'built.pdb',
                                  meshmaker = meshmaker_binary,
                                  radii_file = radii_with_jay_protonation,
                                  resolution = 6)
    # end BEM solver specific stuff

    for titration_state in all_titration_states[]
    
    # restore previous dir
    os.chdir(cwd)
    sys.exit(0)
