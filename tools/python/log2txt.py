#!/usr/bin/env python

# Script:  log2txt.py
# Purpose: extract thermo info from LAMMPS log file
#          create a text file of numbers in columns, suitable for plotting
# Syntax:  log2txt.py log.lammps data.txt X Y ...
#          log.lammps = LAMMPS log file
#          data.txt = text file to create
#          X Y ... = columns to include (optional), X,Y are thermo keywords
#                    if no columns listed, all columns are included
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

from __future__ import print_function

import sys,os,argparse
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.insert(1,path)
from log import log

# set up arg parser
parser = argparse.ArgumentParser()
parser.add_argument('lammpslog', help='name of the lammps log file')
parser.add_argument('outname', help='name of the file to be written')
parser.add_argument('cols', nargs='*', help='any number of column names, optional')
parser.add_argument('-n', action='store_true', help='save column names as the header of the file')

args = parser.parse_args()
logfile = args.lammpslog
datafile = args.outname
columns = args.cols
writenames = args.n

lg = log(logfile)
if columns == []:
  lg.write(datafile, writenames)
else:
  str = "lg.write(datafile, %r" % writenames
  for word in columns: str += ',"' + word + '"'
  str += ')'
  eval(str)
