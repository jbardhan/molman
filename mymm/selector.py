import os, sys, subprocess, argparse, re, logging

import mymm

class Selector:
    def __init__(self, string = None, selector = None, segids = None, resnums = None):
        if string is not None:
            self.string = string.copy()
            return
        elif selector is not None:
            self.string = selector.string.copy()
            return

        if segids is None:
            segids_string =  ''
        else:
            segids_processed = ['(re.match(\'%s\',x.segid) is not None)'%(segid) for segid in segids]
            segids_string = '(' + ' or '.join(segids_processed) + ')'

        string_list = [segids_string]

        if resnums is None:
            resnums_string = ''
        else:
            resnums_string = '(' + 'x.resnum in %s '%(str(resnums))+ ')'
            string_list.append(resnums_string)
            
        self.string = ' and '.join(string_list)

    def __str__(self):
        return selector.string.copy()
