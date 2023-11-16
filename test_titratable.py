#!/usr/bin/env python
import re, os, sys, subprocess, argparse, logging
import mymm

args = {'def_file': 'titratable.def'}
table = mymm.Titratable(filename = args['def_file'])
table.print_titratable_definitions()
