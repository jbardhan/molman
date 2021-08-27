#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise

parser = argparse.ArgumentParser(description="This program takes .pqr and .sites files, and generates the background charge vector and individual charged sites vectors.", prog = sys.argv[0])
parser.add_argument('--srf', metavar = 'lysozyme_4.srf')
parser.add_argument('--dot', metavar = "lysozyme.dot")
parser.add_argument('--nocav', action = 'store_true', default = False)
parser.add_argument('--nostern', action = 'store_true', default = True)
args = parser.parse_args(sys.argv[1:])

srf = mymm.AltmanSrf(args.srf)

