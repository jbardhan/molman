#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

parser = argparse.ArgumentParser(description="This program takes an ICE-style .crd and run.file and converts them into a set of .crg and .crd files corresponding to the desired calculations in the run.file.", prog = sys.argv[0])
parser.add_argument('--crdIn', metavar = 'lysozyme.crd')
parser.add_argument('--runFile', metavar = "run.file")
parser.add_argument('--baseOut', metavar = "lysozyme_", default = "")
args = parser.parse_args(sys.argv[1:])

system = mymm.Molecule()
system.read_crd(args.crdIn)
system.write_crd("checkfile")

runfile = mymm.RunFile(filename = args.runFile, base = args.baseOut)
runfile.print_details()


for curCalc in list(runfile.calculations.keys()):
    modifiedSystem = runfile.process(curCalc, system)
    modifiedSystem.write_crg(runfile.base_output(curCalc) + ".crg")
    modifiedSystem.write_crd(runfile.base_output(curCalc) + ".crd")
    
