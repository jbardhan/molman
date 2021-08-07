import numpy as np
import sys

class Atom:
    def __init__(self, number, atomid, resid, segid, resnum, absres, x, y, z, q=0.0, elem = "", tempfactor = 0.0,occupancy=0.0,radius=0.0):
        self.number = number
        self.atomid = atomid
        self.resid = resid
        self.segid = segid
        self.resnum = resnum
        self.absres = absres
        self.xyz = [float(x), float(y), float(z)]
        self.elem = elem
        self.q   = float(q)
        self.radius     = 0.0
        self.occupancy  = occupancy
        self.tempfactor = tempfactor
        self.radius = radius
        
    def change_resid(self, newresid):
        self.resid = newresid

    def set_number(self, newnumber):
        self.number = newnumber

    def set_radius(self, radii, patches):
        successfully_set_radius = False

        # check if resid is in radii directly
        if self.resid in radii.data:
            if self.atomid in radii.data[self.resid]:
                self.radius = radii.data[self.resid][self.atomid]
                successfully_set_radius = True
        
        # note that in the above clause, we don't just return because we might have a patch that adjusts the radius!!
        print "segid = " + self.segid + " and resnum = " + str(self.resnum)
        patch = patches.get_patch(segid = self.segid, resnum = self.resnum)
        if patch is not None:
            print "Checking patch " + self.segid + " : " + str(self.resnum)
            if self.atomid in radii.data["global"]:
                self.radius = radii.data["global"][self.atomid]
                successfully_set_radius = True

#        if successfully_set_radius == False:
#            print 'Error in atom.set_radius!  Atom atomnum does not have a radius entry anywhere!\n'
#            sys.exit(1)
        
    def set_charge(self, q):
        self.q = q

    def set_q_from_psf(self, psf):
        for atom in psf['ATOMS']:
#            print 'in set_q_from_psf with'
#            print '\tatomid='+str(atom.atomid)+' resnum='+str(atom.resnum)+' resid='+str(atom.resid)+'  ==> q='+str(atom.q)
#            print '\tatomid='+str(self.atomid)+' resnum='+str(self.resnum)+' resid='+str(self.resid)
            
            if self.atomid == atom.atomid and self.resnum == atom.resnum and self.resid == atom.resid:
#                print 'atomid='+str(atom.atomid)+' resnum='+str(atom.resnum)+' resid='+str(atom.resid)+'  ==> q='+str(atom.q)
                self.q = atom.q

    def set_q_from_topology(self, topology, patches):
        successfully_set_charge = False

        if self.resid in topology.residues:
            if self.atomid in topology.residues[self.resid]:
                self.q = topology.residues[self.resid][self.atomid]['charge']
                successfully_set_charge = True

        # Note that we don't "return" from the above clause because we might have a patch that we need to apply!
            
        # ... check if this atom is part of a patch
        patch = patches.get_patch(segid = self.segid, resnum = self.resnum)
        if patch in topology.patches:
            if self.atomid in topology.patches[patch]:
                self.q = topology.patches[patch][self.atomid]['charge']
                successfully_set_charge = True

#        if successfully_set_charge == False:
#            print 'Error in atom.set_q!  Atom atomnum does not have a charge entry anywhere!\n'
#            sys.exit(1)
            
    def crg_entry(self):
        entry_string = "%-6s%3s%4d%1s%8.5f\n"%(self.atomid,self.resid,self.resnum,self.segid,self.q)
        return entry_string
    
    def crd_entry(self):
        entry_string = "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4.4s %-4.4s%10.5f\n"%(self.number, self.absres,self.resid,self.atomid, self.xyz[0], self.xyz[1], self.xyz[2], self.segid, self.resnum, self.q)
        return entry_string
    
    def pdb_entry(self):
        entry_string = "ATOM  %5d %-4.4s %3.3s  %4d    %8.3f%8.3f%8.3f  1.00%6.2f\n" % (self.number, self.atomid, self.resid, self.resnum, self.xyz[0], self.xyz[1], self.xyz[2], self.q)
          
        return entry_string

    # needs better checking against Atom.pm (line 390)
    def pdb2_entry(self):
        atom_string = self.atomid
        if len(atom_string) != 4:
            atom_string = " " + atom_string

        formalcharge = 0
            
        entry_string = "ATOM  %5d %-4.4s %-4.4s%1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4.4s%+2.2s%+2.2s\n" % (self.number, atom_string, self.resid, self.segid[0], self.resnum,self.xyz[0],self.xyz[1],self.xyz[2],self.occupancy,self.tempfactor, self.segid,self.elem,formalcharge)

        return entry_string

    def mead_pqr_entry(self):
        atom_string = self.atomid
        if len(atom_string) != 4:
            atom_string = " " + atom_string

        entry_string = "ATOM  %5d %-4.4s %-4.4s %4d    %8.3f%8.3f%8.3f %7.4f %7.4f\n" % (self.number, atom_string, self.resid, self.resnum,self.xyz[0],self.xyz[1],self.xyz[2],self.q,self.radius)
        return entry_string
    
    def apbs_pqr_entry(self):
        atom_string = self.atomid
        if len(atom_string) != 4:
            atom_string = " " + atom_string

        entry_string = "ATOM  %5d %-4.4s %-4.4s %-2.2s %4d    %8.3f%8.3f%8.3f %6.3f %6.3f\n" % (self.number, atom_string, self.resid, self.segid, self.resnum,self.xyz[0],self.xyz[1],self.xyz[2],self.q,self.radius)
        return entry_string

    def accurate_pqr_entry(self):
        atom_string = self.atomid
        if len(atom_string) != 4:
            atom_string = " " + atom_string

        entry_string = "ATOM  %5d %-4.4s %-4.4s %-2.2s %4d    %10.6f%10.6f%10.6f %10.6f %10.6f\n" % (self.number, atom_string, self.resid, self.segid, self.resnum,self.xyz[0],self.xyz[1],self.xyz[2],self.q,self.radius)
        return entry_string


    def xyzr_entry(self):
        entry_string = "%8.3f %8.3f %8.3f %6.3f\n" %(self.xyz[0],self.xyz[1],self.xyz[2],self.radius)
        return entry_string

    def rotateZ(self, angle = 0.0):
        # note: angle is in degrees, not radians
        angle = angle * np.pi / 180.0
        sinTheta = np.sin(angle)
        cosTheta = np.cos(angle)
        matrix  = np.array([
            [ cosTheta, -sinTheta, 0.0],
            [ sinTheta, cosTheta, 0.0],
            [ 0.0, 0.0, 1.0]])
        xyz = np.array([
            [self.xyz[0]],
            [self.xyz[1]],
            [self.xyz[2]]])
        newxyz = np.dot(matrix,xyz)
#        print("Newxyz is " + str(newxyz))
        self.xyz[0] = float(newxyz[0])
        self.xyz[1] = float(newxyz[1])
        self.xyz[2] = float(newxyz[2])
        
    def calculate_dipole_moment(self):
        return [self.xyz[0]*self.q,self.xyz[1]*self.q,self.xyz[2]*self.q]
    
    def translate(self, vec = None, x = None, y = None, z = None):
        if (vec == None) and (x==None or y==None or z==None):
            print "Error: translate needs either a vector or all x,y,z!\n"
            sys.exit(1)

        if (vec != None and x != None and y != None and z != None):
            print "Error: translate needs _either_ vec or x,y,z, but not both!\n"
            sys.exit(1)
            
        if vec != None:
            x = vec[0]
            y = vec[1]
            z = vec[2]
            
        self.xyz[0] = self.xyz[0] + x
        self.xyz[1] = self.xyz[1] + y
        self.xyz[2] = self.xyz[2] + z
