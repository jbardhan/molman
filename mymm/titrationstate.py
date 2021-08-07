class TitrationState:
    def __init__(self, filename):
        self.segids_titrating = {}
        fh = open(filename,'r')
        for line in fh:
            line_data = line.rstrip().lstrip().split()

            if line_data[0] not in self.segids_titrating.keys():
                self.segids_titrating[line_data[0]] = {}

            if line_data[1] not in self.segids_titrating[line_data[0]].keys():
                self.segids_titrating[line_data[0]][line_data[1]] = line_data[2]

    def atom_is_in_titrating_group(self, mutator, atom):
        if atom.segid not in self.segids_titrating.keys():
            return 0

        if str(atom.resnum) not in self.segids_titrating[atom.segid].keys():
            return 0
        
        return 1

    def print_state(self):
        print "TitrationState object:"
        print "SegIDs with titrating groups: " + " ".join(self.segids_titrating.keys())
        for cursegid in self.segids_titrating.keys():
            print "SegID " + cursegid + ": " + " ".join(self.segids_titrating[cursegid].keys())

