#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm


args = { 'pdb': 'arg_capped_nter_cter_PARSE.pdb',
         'psf': 'arg_capped_nter_cter_PARSE.psf',
         'def_file': '../titratable_PARSE.def',
         'titration_list' : 'titration_list.txt',
         'topology' : '/home/bard415/repos/60376-testset1/top_PARSE_prot_patch.inp',
         'radii' : '/home/bard415/repos/parameters/radii_fixNH3.siz'
}


protein = mymm.Apbs()
protein.pqrList.append("neutral_protein.pqr")
protein.headerComment = "Junk test APBS input"


protein_params_hash =  {'dime': 33,
                        'cglen': [40, 40, 40],
                        'fglen': [20, 20, 20],
                        'fgcent': 1,
                        'cgcent': 1}

protein.addElecSection(name="solv", molIndex=1, commands = "write atompot flat solv", additionalParamsDict=protein_params_hash)

protein_params_hash.update({'sdie':protein.elecParams['pdie']})

protein.addElecSection(name="ref", molIndex=1, commands = "write atompot flat ref", additionalParamsDict=protein_params_hash)

protein.addAnalysisSection(["print elecEnergy solv - ref end"])
protein.printApbsInput("foo.in")
