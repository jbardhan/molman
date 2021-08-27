import ntpath, re, sys

class AltmanSrf:
    def __init__(self, filename):
        self.filename = filename;
        self.srf_path, self.srf_file = ntpath.split(filename)

        fh = open(filename,'r')
        stern_panel_type = fh.readline().rstrip().lstrip()
        if stern_panel_type == 'c':
            self.stern_panel_type = 'curved'
        elif stern_panel_type == 'f':
            self.stern_panel_type = 'flat'
        else:
            print(("Error in AltmanSrf: file " + filename + " has " + stern_panel_type + " as a panel type for\n"))
            print("\tthe Stern layer surfaces.  Only 'c' and 'f' are allowed.  Dying...\n")
            sys.exit(0)
        print(("Stern panel type is " + self.stern_panel_type))
        
        dielectric_panel_type = fh.readline().rstrip().lstrip()
        if dielectric_panel_type == 'c':
            self.dielectric_panel_type = 'curved'
        elif dielectric_panel_type == 'f':
            self.dielectric_panel_type = 'flat'
        else:
            print(("Error in AltmanSrf: file " + filename + " has " + dielectric_panel_type + " as a panel type for\n"))
            print("\tthe outer dielectric surfaces.  Only 'c' and 'f' are allowed.  Dying...\n")
            sys.exit(0)
        print(("Dielectric panel type is " + self.dielectric_panel_type))

        outer_stern_file_root = fh.readline().rstrip().lstrip()
        self.outer_stern_file_root = self.srf_path + "/./" + outer_stern_file_root

        print(("outer Stern layer file root, full path, is " + self.outer_stern_file_root))
        
        outer_dielectric_file_root = fh.readline().rstrip().lstrip()
        self.outer_dielectric_file_root = self.srf_path + "/./" + outer_dielectric_file_root

        dielectric_cavity_file_root = fh.readline().rstrip().lstrip()
        stern_cavity_file_root = fh.readline().rstrip().lstrip()
        outer_dielectric_inside = fh.readline().rstrip().lstrip().split()
        cavity_inside = fh.readline().rstrip().lstrip().split()
        stern_cavity_inside = fh.readline().rstrip().lstrip().split()
        fh.close()
        
    
