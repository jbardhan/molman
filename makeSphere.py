#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging
import mymm

sph = mymm.spheresolute.create_sphere_of_atoms(radius = 7.2, density = 1.6)
sph.write_pdb2("sph8.pdb")

mymol = mymm.Molecule()
mymol.add_atom(mymm.Atom(number=1, atomid='CA',resid='SPH',segid='F',resnum=991,absres="",x = 0.0,y=0.0,z=0.0,q=1.))
mymol.add_atom(mymm.Atom(number=2, atomid='CA',resid='SPH',segid='A',resnum=992,absres="",x = 0.0,y=0.0,z=7.0,q=1.0))
mymol.add_atom(mymm.Atom(number=3, atomid='CA',resid='SPH',segid='A',resnum=992,absres="",x = 0.0,y=0.0,z=5.0,q=-1.0))
mymol.add_atom(mymm.Atom(number=4, atomid='CA',resid='SPH',segid='B',resnum=993,absres="",x = 7.0,y=0.0,z=0.0,q=-1.0))
mymol.add_atom(mymm.Atom(number=5, atomid='CA',resid='SPH',segid='B',resnum=993,absres="",x = 5.0,y=0.0,z=0.0,q=1.0))

mymol.write_pdb2("molecule.pdb")

              
