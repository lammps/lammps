#!/usr/bin/env python -i

# Script:  logplot.py
# Purpose: use GnuPlot to plot two columns from a LAMMPS log file
# Syntax:  logplot.py log.lammps X Y
#          log.lammps = LAMMPS log file
#          X,Y = plot Y versus X where X,Y are thermo keywords
#          once plot appears, you are in Python interpreter, type C-D to exit
# Author:  Steve Plimpton (Sandia), sjplimp at sandia.gov

import sys,os
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.append(path)
from log import log
from gnu import gnu

if len(sys.argv) != 4:
  raise StandardError, "Syntax: logplot.py log.lammps X Y"

logfile = sys.argv[1]
xlabel = sys.argv[2]
ylabel = sys.argv[3]

lg = log(logfile)
x,y = lg.get(xlabel,ylabel)
g = gnu()
g.plot(x,y)
print "Type Ctrl-D to exit Python"
