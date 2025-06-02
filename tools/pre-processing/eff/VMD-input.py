#!/usr/local/bin/python-2.5/bin/python

import sys, os
from getopt import gnu_getopt as getopt

Info="""
Module name: VMD-input.py 

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@wag.caltech.edu
Project: pEFF
Version: August 2009

Usage: python VMD-input.py lammps_dump_filename radii_column xpos_column
Example: python VMD-input.py dump.lammpstrj 5 6

1. Extracts the electron radii from a lammps trajectory dump into %s.out
2. Creates %s.xyz file
3. Creates appropriate %.vmd file which can be sourced using TCL/TK in VMD

"""

from lmp2radii_col import makeradii 
from lmp2xyz import lmp2xyz
 
def printHelp(input):
  Info%(input,input,input)

if __name__ == '__main__':
    
   # if no input, print help and exit
    if len(sys.argv) < 2:
        print "Usage: python VMD-input.py lammps_dump_filename radii_column xpos_column\n"
        sys.exit(1)
    else:
        infile=sys.argv[1]

    workdir=os.getcwd()
    tools=sys.argv[0].split('VMD-input.py')[0]
    # set defaults
    outfile = infile.split('.')[0]
    print sys.argv
    if len(sys.argv) == 4:
      column = int(sys.argv[2])
      xpos = int(sys.argv[3])
      print "Assuming xpos=eradius+1"
    elif len(sys.argv) == 3:
      column = int(sys.argv[2])
      xpos = column+1
      print "Assuming xpos=eradius+1"
    elif len(sys.argv) == 2:
      column=5    # default = radius for dump -> id type q spin eradius x y z
      xpos=6
    else: print "Incorrect number of arguments"

    # check for input:
    opts, argv = getopt(sys.argv[1:], 'c:o:ha')

    # read options
    for opt, arg in opts:
        if opt == '-h':             # -h: print help
          Info%(input,input,input)
        if opt == '-o':         # output file name
          outfile=arg
        if opt == '-c':         # select column from lammpstrj file to tabulate
          column=int(arg)
    print column,xpos
    makeradii(infile,outfile+".out",column,True)
    lmp2xyz(infile,outfile+".xyz",xpos)
    print "Creating %s file ..."%(outfile+".vmd")
    os.system("cat %s | sed 's/xyzfile/%s/' > %s"%(tools+"radii.vmd",outfile+".xyz","temp"))
    os.system("cat %s | sed 's/radiifile/%s/' > %s; rm temp"%("temp",outfile+".out",workdir+'/'+outfile+".vmd"))
    print "Done !! (you can now source %s using VMD's console) \n"%(outfile+".vmd")

    print "NOTE: In VMD, set graphics representation for electrons to transparency,"
    print "and change the atom types in the xyz file according to your values,"
    print "for simplicity, they are set using the same mass sequence definition\nfrom your lammps data file\n"
