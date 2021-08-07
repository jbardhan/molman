class TitrationData:
    def __init__(self, filename):
        self.residue_data = {}
        self.read_residue_data(filename)

    def read_residue_data(self,filename):
        titration_file = open(args.titr)

        for line in titration_file:
            line_comment = re.search("#", line)
            if line_comment is not None:
                line = line[0:line_comment.start()]

            line_init_residue = re.search("TITR", line)
            if line_init_residue is not None:
                line_init_data = line.rstrip().lstrip().split()
                self.residue_data[line_init_data[1]] = {'atoms':{},'pKa':line_init_data[2],'numProtonationStates':line_init_data[3],'protons':line_init_data[4:]}
                continue
        
            line_data = line.rstrip().lstrip().split()

            if len(line_data) == 0:
                continue
    
            if line_data[0] not in self.residue_data.keys():
                print "Error: You must define the residue first with TITR <resname> <pka> <protons>\n"
                sys.exit(2)
        #        self.residue_data[line_data[0]] = {'atoms':{},'pKa':None,'numProtonationStates':0,'protons':[]}

            if line_data[1] in self.residue_data[line_data[0]]['atoms'].keys():
                print "Duplicate atom " + line_data[1] + " in residue " + line_data[0]
                sys.exit(2)
        
            last_prot_charge_index = len(line_data)-1;
            self.residue_data[line_data[0]]['atoms'][line_data[1]] = {'charges':line_data[3:last_prot_charge_index]}

        titration_file.close()

