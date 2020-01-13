# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# cfg tool

oneline = "Convert LAMMPS snapshots to AtomEye CFG format"

docstr = """
c = cfg(d)              d = object containing atom coords (dump, data)

c.one()                 write all snapshots to tmp.cfg
c.one("new")            write all snapshots to new.cfg
c.many()                write snapshots to tmp0000.cfg, tmp0001.cfg, etc
c.many("new")           write snapshots to new0000.cfg, new0001.cfg, etc
c.single(N)             write snapshot for timestep N to tmp.cfg
c.single(N,"file")      write snapshot for timestep N to file.cfg
"""

# History
#   11/06, Aidan Thompson (SNL): original version

# ToDo list
# should decide if dump is scaled or not, since CFG prints in scaled coords
# this creates a simple AtomEye CFG format
#   there is more complex format we could write out
#   which allows for extra atom info, e.g. to do atom coloring on
# how to dump for a triclinic box, since AtomEye accepts this

# Variables
#   data = data file to read from

# Imports and external programs

import sys

# Class definition

class cfg:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data

  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) == 0: file = "tmp.cfg"
    elif args[0][-4:] == ".cfg": file = args[0]
    else: file = args[0] + ".cfg"

    f = open(file,"w")
    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      xlen = box[3]-box[0]
      ylen = box[4]-box[1]
      zlen = box[5]-box[2]

      print >>f,"Number of particles = %d " % len(atoms)
      print >>f,"# Timestep %d \n#\nA = 1.0 Angstrom" % time
      print >>f,"H0(1,1) = %20.10f A " % xlen
      print >>f,"H0(1,2) = 0.0 A "
      print >>f,"H0(1,3) = 0.0 A "
      print >>f,"H0(2,1) = 0.0 A "
      print >>f,"H0(2,2) = %20.10f A " % ylen
      print >>f,"H0(2,3) = 0.0 A "
      print >>f,"H0(3,1) = 0.0 A "
      print >>f,"H0(3,2) = 0.0 A "
      print >>f,"H0(3,3) = %20.10f A " % zlen
      print >>f,"#"

      for atom in atoms:
        itype = int(atom[1])
        xfrac = (atom[2]-box[0])/xlen
        yfrac = (atom[3]-box[1])/ylen
        zfrac = (atom[4]-box[2])/zlen
#        print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f " % (itype,xfrac,yfrac,zfrac,atom[5],atom[6],atom[7])
        print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  0.0 0.0 0.0 " % (itype,xfrac,yfrac,zfrac)

      print time,
      sys.stdout.flush()
      n += 1

    f.close()
    print "\nwrote %d snapshots to %s in CFG format" % (n,file)

  # --------------------------------------------------------------------

  def many(self,*args):
    if len(args) == 0: root = "tmp"
    else: root = args[0]

    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      if n < 10:
        file = root + "000" + str(n)
      elif n < 100:
        file = root + "00" + str(n)
      elif n < 1000:
        file = root + "0" + str(n)
      else:
        file = root + str(n)
      file += ".cfg"
      f = open(file,"w")

      xlen = box[3]-box[0]
      ylen = box[4]-box[1]
      zlen = box[5]-box[2]

      print >>f,"Number of particles = %d " % len(atoms)
      print >>f,"# Timestep %d \n#\nA = 1.0 Angstrom" % time
      print >>f,"H0(1,1) = %20.10f A " % xlen
      print >>f,"H0(1,2) = 0.0 A "
      print >>f,"H0(1,3) = 0.0 A "
      print >>f,"H0(2,1) = 0.0 A "
      print >>f,"H0(2,2) = %20.10f A " % ylen
      print >>f,"H0(2,3) = 0.0 A "
      print >>f,"H0(3,1) = 0.0 A "
      print >>f,"H0(3,2) = 0.0 A "
      print >>f,"H0(3,3) = %20.10f A " % zlen
      print >>f,"#"

      for atom in atoms:
        itype = int(atom[1])
        xfrac = (atom[2]-box[0])/xlen
        yfrac = (atom[3]-box[1])/ylen
        zfrac = (atom[4]-box[2])/zlen
#        print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f " % (itype,xfrac,yfrac,zfrac,atom[5],atom[6],atom[7])
        print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  0.0 0.0 0.0 " % (itype,xfrac,yfrac,zfrac)

      print time,
      sys.stdout.flush()
      f.close()
      n += 1

    print "\nwrote %s snapshots in CFG format" % n

  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) == 0: file = "tmp.cfg"
    elif args[0][-4:] == ".cfg": file = args[0]
    else: file = args[0] + ".cfg"

    which = self.data.findtime(time)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    f = open(file,"w")

    xlen = box[3]-box[0]
    ylen = box[4]-box[1]
    zlen = box[5]-box[2]

    print >>f,"Number of particles = %d " % len(atoms)
    print >>f,"# Timestep %d \n#\nA = 1.0 Angstrom" % time
    print >>f,"H0(1,1) = %20.10f A " % xlen
    print >>f,"H0(1,2) = 0.0 A "
    print >>f,"H0(1,3) = 0.0 A "
    print >>f,"H0(2,1) = 0.0 A "
    print >>f,"H0(2,2) = %20.10f A " % ylen
    print >>f,"H0(2,3) = 0.0 A "
    print >>f,"H0(3,1) = 0.0 A "
    print >>f,"H0(3,2) = 0.0 A "
    print >>f,"H0(3,3) = %20.10f A " % zlen
    print >>f,"#"

    for atom in atoms:
      itype = int(atom[1])
      xfrac = (atom[2]-box[0])/xlen
      yfrac = (atom[3]-box[1])/ylen
      zfrac = (atom[4]-box[2])/zlen
#        print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f " % (itype,xfrac,yfrac,zfrac,atom[5],atom[6],atom[7])
      print >>f,"1.0  %d   %15.10f  %15.10f  %15.10f  0.0 0.0 0.0 " % (itype,xfrac,yfrac,zfrac)

    f.close()
