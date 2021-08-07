import os, sys, subprocess, argparse, re, logging
import math
import numpy as np
import mymm


def get_random_points_on_sphere(radius, numpoints):
    t = np.arange(start=-1.0+1.0/numpoints, stop=1.0+1.0/numpoints, step=2.0/numpoints)
    theta = np.arccos(t)
    phi = math.sqrt(math.pi*numpoints)*np.arcsin(t)
    points = radius* np.array([np.sin(theta)*np.cos(phi),
                               np.sin(theta)*np.sin(phi),
                               np.cos(theta)])
    return list(points.transpose())

def create_sphere_of_atoms(radius, density = 0.25, origin = [0.0, 0.0, 0.0], segid = "S"):
    newsphere = mymm.Molecule()
    if radius < 2:
        print "Error: addSphereOfAtoms only works for radius > 2!\n"
        return newsphere

    atomid = "CA"
    resid  = "SPH"
    number = 1
    resnum = 1
    radius_for_atom_set = radius-2.0
    while radius_for_atom_set >= 2.0:
        surface_area = 4.0 * math.pi * radius_for_atom_set**2
        num_points = int(surface_area * density)
        print "radius = %f, num_points =  %d\n" % (radius_for_atom_set, num_points)
        points = get_random_points_on_sphere(radius_for_atom_set,num_points)
        for curpoint in points:
            newsphere.add_atom(mymm.Atom(number, atomid, resid, segid,
                                   resnum, resnum,
                                   curpoint[0]+origin[0],
                                   curpoint[1]+origin[1],
                                   curpoint[2]+origin[2]))
            number = number + 1
            resnum = resnum + 1
        radius_for_atom_set = radius_for_atom_set - 3.0

    return newsphere
