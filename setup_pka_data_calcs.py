#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm
global_param_dir = os.getenv('MOLMAN_PARAMETER_DIR',os.getcwd())
#proj_param_dir   = os.getenv('HOME') + "/repos/parameters"

args = {'pdb':'built.pdb',
        'psf':'built.psf',
        'pqr':'built.pqr',
        'group_defs':'titratable_PARSE.def',
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

table = mymm.Titratable(filename = args['group_defs'])

titratable.read_titration_list(args['titratable_list'])

print "Titratable residues are: " + str(titratable.list_of_residues)

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
