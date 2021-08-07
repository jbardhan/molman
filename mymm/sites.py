import sys

class Sites:
    def __init__(self, filename):
        self.site_data = []
        self.site_types = {}
        
        fh = open(filename, 'r')
        for line in fh:
            line_data = line.rstrip().lstrip().split()
            if len(line_data) != 2:
                print "Error!  Sites file had a line with more than two fields: " + line
                sys.exit(1)
            (resnum, site_type) = tuple(line_data)
            self.site_data.append({'resnum':int(resnum), 'type':site_type})
            if site_type not in self.site_types:
                self.site_types[site_type] = {}
        fh.close()

        for site_type in self.site_types.keys():
            site_filename = site_type + ".st"
            fh = open(site_filename,'r')

            atom_data = []
            pka_line = fh.readline()
            pka_model = float(pka_line.rstrip().lstrip())
            total_protonated_charge = 0.0
            total_unprotonated_charge = 0.0

            for line in fh:
                line_data = line.rstrip().lstrip().split()
                if len(line_data) != 4:
                    print "Error! site file " + site_filename + " has a line with more than 4 fields...\n\t" + line
                    sys.exit(1)
                (resid, atomid, protonated_q, unprotonated_q) = tuple(line_data)
                protonated_q = float(protonated_q)
                unprotonated_q = float(unprotonated_q)
                total_protonated_charge   += protonated_q
                total_unprotonated_charge += unprotonated_q
                atom_data.append({'atomid':atomid,'prot':protonated_q,'unprot':unprotonated_q})
            fh.close()

            if abs(total_unprotonated_charge - int(total_unprotonated_charge)) > 1e-3:
                print "Error! Total unprotonated charge on site " + site_type + ", resid " + resid + ", is " + str(total_unprotonated_charge) + "."
                sys.exit(1)
            if abs(total_protonated_charge - int(total_protonated_charge)) > 1e-3:
                print "Error! Total protonated charge on site " + site_type + ", resid " + resid + ", is " + str(total_protonated_charge) + "."
                sys.exit(1)
            if (total_protonated_charge < total_unprotonated_charge) or abs(total_protonated_charge - total_unprotonated_charge - 1.0) > 1e-3:
                print "Error! Total protonated charge " + str(total_protonated_charge) + " is not 1.0 greater than total unprotonated charge " + str(total_unprotonated_charge)
                sys.exit(1)

            if abs(total_protonated_charge - 1) < 1e-3:
                print "Site " + site_type + ", resid " + resid + " is assumed to be basic, since the total protonated charge is " + str(total_protonated_charge)
                acidic = False
            else:
                print "Site " + site_type + ", resid " + resid + " is assumed to be acidic, since the total protonated charge is " + str(total_protonated_charge)
                acidic = True
                
            self.site_types[site_type] = {'resid':resid, 'pKa_model':pka_model, 'atoms':atom_data, 'acidic':acidic}

    def set_zero_charges(self, molecule):
        for atom in molecule.atoms:
            atom.set_charge(0.0)
            
    def set_protonated(self, molecule, site):
        if self.is_acidic(site):
            self.set_uncharged(molecule, site)
        else:
            self.set_charged(molecule,site)

    def set_deprotonated(self, molecule, site):
        if self.is_acidic(site):
            self.set_charged(molecule, site)
        else:
            self.set_uncharged(molecule,site)
            
    def set_neutral(self, molecule):
        for cur_site in self.site_data:
            self.set_uncharged(molecule, cur_site)

    def set_uncharged(self, molecule, site):
        if self.is_acidic(site):
            self.set_charges_on_site(molecule, site, 'protonated')
        else:
            self.set_charges_on_site(molecule, site, 'unprotonated')

    def set_charged(self, molecule, site):
        if self.is_acidic(site):
            self.set_charges_on_site(molecule, site, 'unprotonated')
        else:
            self.set_charges_on_site(molecule, site, 'protonated')

    def is_acidic(self, site):
        return self.site_types[site['type']]['acidic']

    def set_charges_on_site(self, molecule, site, charge_state):
#        print "Setting charges on resnum " + str(site['resnum']) 
        
        for atom_to_match in self.site_types[site['type']]['atoms']:
            atom = molecule.find_atom(resnum = site['resnum'], resid = self.site_types[site['type']]['resid'], atomid = atom_to_match['atomid'])

            if charge_state == 'protonated':
                atom.set_charge(atom_to_match['prot'])
#                print "setting charge on " + atom.atomid + " to protonated"
            elif charge_state == 'unprotonated':
                atom.set_charge(atom_to_match['unprot'])
#                print "setting charge to unprotonated"
            else:
                print "Error!  charge_state " + charge_state + " is not defined!"
                sys.exit(1)
            
           
    
