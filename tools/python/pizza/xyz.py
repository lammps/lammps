# Pizza.py toolkit, www.cs.sandia.gov/~sjplimp/pizza.html
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# xyz tool

oneline = "Convert LAMMPS snapshots to XYZ format"

docstr = """
x = xyz(d)              d = object containing atom coords (dump, data)

x.one()                 write all snapshots to tmp.xyz
x.one("new")            write all snapshots to new.xyz
x.many()                write snapshots to tmp0000.xyz, tmp0001.xyz, etc
x.many("new")           write snapshots to new0000.xyz, new0001.xyz, etc
x.single(N)             write snapshot for timestep N to tmp.xyz
x.single(N,"file")      write snapshot for timestep N to file.xyz
"""

# History
#   8/05, Steve Plimpton (SNL): original version

# ToDo list

# Variables
#   data = data file to read from

# Imports and external programs

import sys

# Class definition

class xyz:

  # --------------------------------------------------------------------

  def __init__(self,data):
    self.data = data

  # --------------------------------------------------------------------

  def one(self,*args):
    if len(args) == 0: file = "tmp.xyz"
    elif args[0][-4:] == ".xyz": file = args[0]
    else: file = args[0] + ".xyz"

    f = open(file,"w")
    n = flag = 0
    while 1:
      which,time,flag = self.data.iterator(flag)
      if flag == -1: break
      time,box,atoms,bonds,tris,lines = self.data.viz(which)

      print >>f,len(atoms)
      print >>f,"Atoms"
      for atom in atoms:
        itype = int(atom[1])
        print >>f,itype,atom[2],atom[3],atom[4]

      print time,
      sys.stdout.flush()
      n += 1

    f.close()
    print "\nwrote %d snapshots to %s in XYZ format" % (n,file)

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
      file += ".xyz"
      f = open(file,"w")
      print >>f,len(atoms)
      print >>f,"Atoms"
      for atom in atoms:
        itype = int(atom[1])
        print >>f,itype,atom[2],atom[3],atom[4]
      print time,
      sys.stdout.flush()
      f.close()
      n += 1

    print "\nwrote %s snapshots in XYZ format" % n

  # --------------------------------------------------------------------

  def single(self,time,*args):
    if len(args) == 0: file = "tmp.xyz"
    elif args[0][-4:] == ".xyz": file = args[0]
    else: file = args[0] + ".xyz"

    which = self.data.findtime(time)
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    f = open(file,"w")
    print >>f,len(atoms)
    print >>f,"Atoms"
    for atom in atoms:
      itype = int(atom[1])
      print >>f,itype,atom[2],atom[3],atom[4]
    f.close()
