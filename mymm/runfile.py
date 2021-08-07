import re
import mymm

def process_line(line):
    import re
    return re.findall(r'\"(.+?)\"',line)

class RunFile:
    def __init__(self, filename, base = ""):
        self.calculations = {}
        self.userBase = base

        fh = open(filename, 'r')
        inCalc = False  # in BEM only have to worry about the ones marked as "final"
        newCalc = False
        for line in fh:
            if (line.find('mark=final')==0) or (line.find('mark=reference')==0):
                newCalc = True

            if newCalc and (line.find('name_group')==0):
                newCalc = False
                inCalc = True
                name_group = process_line(line)[0]  # process_line returns a list
                self.calculations[name_group] = {}

            if inCalc and (line.find('atoms_charged')==0):
                self.calculations[name_group]['charge_list'] = process_line(line)

            if inCalc and (line.find('atoms_shape')==0):
                self.calculations[name_group]['shape_list'] = process_line(line)
                
            
        fh.close()

    def process(self, current, system):
        modified_system = mymm.Molecule()
        for atom in system.atoms:
            if atom.segid in self.calculations[current]['shape_list']:
                modified_system.add_atom(atom)
                
        for atom in modified_system.atoms:
            if atom.segid not in self.calculations[current]['charge_list']:
                atom.set_charge(q = 0.0)

        return modified_system

    def base_output(self, current):
        return self.userBase + current
            
    def print_details(self):
        print "RunFile:"
        print "\tuserBase = " + self.userBase

        for key,val in self.calculations.iteritems():
            print "\tname_group = " + key
            for key2,val2 in val.iteritems():
                print "\t\t" + key2 + ": " + " ".join(val2)

    
