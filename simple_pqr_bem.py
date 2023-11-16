#!/usr/bin/env python

import os, sys, subprocess, argparse, re, logging, errno
import mymm
import tempfile, shutil
import ntpath 


def write_delphi_param_file(outputfilename, details):
    f = open(outputfilename, "w")
    f.write("DVERT="+str(details.dieldens)+"\n")
    f.write("SVERT="+str(details.sterndens)+"\n")
    f.write("INDI="+str(details.idiel)+"\n")
    f.write("EXDI="+str(details.odiel)+"\n")
    f.write("SALT="+str(details.saltconc)+"\n")
    f.write("IONRAD="+str(details.stern)+"\n")
    f.write("PRBRAD="+str(details.proberad)+"\n")
    f.write("DE=1e-4\n")
    f.write("ENERGY(G)\n")
    f.close()


def mesh_surface(meshmaker_bin, outputXyzr, radiiSiz, outputSrf, proberadius, sternradius, dieldens, sterndens, saltconc, prefix):
#    Usage: /home/bard415/repos/fftsvd/meshmaker2 [molecule.pdb/crd/xyzr] [radii.siz] [output.srf] [proberadius] [ionexclusionradius] [dielectricdensity] [saltdensity] [dielectriccode] [saltcode] [prefix]
    meshmaker_command = ' '.join([meshmaker_bin, outputXyzr, radiiSiz, outputSrf, str(proberadius), str(sternradius), str(dieldens), str(sterndens), str(saltconc), "1", "1", prefix])
    print("Meshmaker command = \n")
    print("\t"+ meshmaker_command + "\n")
#    os.system(meshmaker_command)
    
             


parser = argparse.ArgumentParser(description = "This program takes a .pqr file (MEAD format only for now, meaning no chain field) and writes a CRG and PDB file from it.", prog = sys.argv[0])

parser.add_argument('--pqr', metavar = 'ligand.pqr')
parser.add_argument('--radii', metavar = 'radii.siz')
parser.add_argument('--idiel', metavar = 'eps_prot', type = float, default = 1.0)
parser.add_argument('--odiel', metavar = 'eps_water', type = float, default = 80.0)
parser.add_argument('--saltconc', metavar = 'salt_concentration', type = float, default = 0.145)
parser.add_argument('--stern', metavar = 'Stern layer thickness', type = float, default = 2.0)
parser.add_argument('--usecavities',metavar = 'Use cavities?', default = True)
parser.add_argument('--proberad',metavar = 'Probe radius', type = float, default = 1.4)
parser.add_argument('--dieldens', metavar = 'Vertex density for dielectrics', type = float, default = 4.0)
parser.add_argument('--sterndens',metavar = 'Vertex density for Sterns', type = float, default = 4.0)
parser.add_argument('--outputdir', metavar = '<directory for outputs>')
parser.add_argument('--outputbase', metavar = '<base for output filename>')
args = parser.parse_args(sys.argv[1:])

args.outputdir = os.path.abspath(args.outputdir)

# make temp dir, change into it
f = tempfile.TemporaryDirectory() # see guide here https://www.tutorialspoint.com/generate-temporary-files-and-directories-using-python
f.name = "deletedir"
startingdir = os.getcwd()

# copy PQR and radii.siz files into temp dir
shutil.copyfile(args.pqr, f.name + "/" + ntpath.basename(args.pqr))
localRadiiSizName = f.name+"/" + ntpath.basename(args.radii)
shutil.copyfile(args.radii, localRadiiSizName)


os.chdir(f.name)
print("Now in temp dir "+f.name)

system = mymm.Molecule()
system.read_pqr(args.pqr)


outputPdb = args.outputbase + ".pdb"
outputCrg = args.outputbase + ".crg"

system.write_pdb2(outputPdb)
system.write_crg(outputCrg)

outputXyzr = args.outputbase + ".xyzr"
system.write_xyzr(outputXyzr)
outputPrm = args.outputbase+".prm"
write_delphi_param_file(outputfilename=outputPrm, details = args)

mesh_surface("~/bin/meshmaker", outputXyzr, ntpath.basename(args.radii), "./output.srf" , args.proberad, args.stern, args.dieldens, args.sterndens, args.saltconc, "./")

#run_bem()
#process_bem_output()

# copy files to outputdir
shutil.copyfile(outputPdb, args.outputdir + "/" + outputPdb)
shutil.copyfile(outputCrg, args.outputdir + "/" + outputCrg)
shutil.copyfile(outputPrm, args.outputdir + "/" + outputPrm)


# clean up, go home 
os.chdir(startingdir)
#f.cleanup() # see guide here https://www.tutorialspoint.com/generate-temporary-files-and-directories-using-python


