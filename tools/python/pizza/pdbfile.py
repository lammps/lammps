# Pizza.py toolkit, https://lammps.github.io/pizza
# LAMMPS development team: developers@lammps.org
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# for python3 compatibility
from __future__ import print_function

# pdb tool

oneline = "Read, write PDB files in combo with LAMMPS snapshots"

docstr = """
p = pdbfile("3CRO")         create pdb object from PDB file or WWW
p = pdbfile("pep1 pep2")    read in multiple PDB files
p = pdbfile("pep*")         can use wildcards
p = pdbfile(d)              read in snapshot data with no PDB file
p = pdbfile("3CRO",d)       read in single PDB file with snapshot data

  string arg contains one or more PDB files
    don't need .pdb suffix except wildcard must expand to file.pdb
    if only one 4-char file specified and it is not found,
      it will be downloaded from http://www.rcsb.org as 3CRO.pdb
  d arg is object with atom coordinates (dump, data)

p.one()                     write all output as one big PDB file to tmp.pdb
p.one("mine")               write to mine.pdb
p.many()                    write one PDB file per snapshot: tmp0000.pdb, ...
p.many("mine")              write as mine0000.pdb, mine0001.pdb, ...
p.single(N)                 write timestamp N as tmp.pdb
p.single(N,"new")           write as new.pdb

  how new PDB files are created depends on constructor inputs:
    if no d: one new PDB file for each file in string arg (just a copy)
    if only d specified: one new PDB file per snapshot in generic format
    if one file in str arg and d: one new PDB file per snapshot
      using input PDB file as template
    multiple input PDB files with a d is not allowed

index,time,flag = p.iterator(0)
index,time,flag = p.iterator(1)

  iterator = loop over number of PDB files
    call first time with arg = 0, thereafter with arg = 1
    N = length = # of snapshots or # of input PDB files
    index = index of snapshot or input PDB file (0 to N-1)
    time = timestep value (time stamp for snapshot, index for multiple PDB)
    flag = -1 when iteration is done, 1 otherwise
  typically call p.single(time) in iterated loop to write out one PDB file
"""

# History
#   8/05, Steve Plimpton (SNL): original version
#   3/17, Richard Berger (Temple U): improve Python 3 compatibility

# ToDo list
#   for generic PDB file (no template) from a LJ unit system,
#     the atoms in PDB file are too close together

# Variables
#   files = list of input PDB files
#   data = data object (ccell,data,dump) to read snapshots from
#   atomlines = dict of ATOM lines in original PDB file
#     key = atom id, value = tuple of (beginning,end) of line

# Imports and external programs

import sys, glob, urllib
PY3 = sys.version_info[0] == 3

if PY3:
    string_types = str,
else:
    string_types = basestring

# Class definition

class pdbfile:

  # --------------------------------------------------------------------

  def __init__(self,*args):
    if len(args) == 1:
      if type(args[0]) is string_types:
        filestr = args[0]
        self.data = None
      else:
        filestr = None
        self.data = args[0]
    elif len(args) == 2:
      filestr = args[0]
      self.data = args[1]
    else: raise Exception("invalid args for pdb()")

    # flist = full list of all PDB input file names
    # append .pdb if needed

    if filestr:
      list = filestr.split()
      flist = []
      for file in list:
        if '*' in file: flist += glob.glob(file)
        else: flist.append(file)
      for i in range(len(flist)):
        if flist[i][-4:] != ".pdb": flist[i] += ".pdb"
      if len(flist) == 0:
        raise Exception("no PDB file specified")
      self.files = flist
    else: self.files = []

    if len(self.files) > 1 and self.data:
      raise Exception("cannot use multiple PDB files with data object")
    if len(self.files) == 0 and not self.data:
      raise Exception("no input PDB file(s)")

    # grab PDB file from http://rcsb.org if not a local file

    if len(self.files) == 1 and len(self.files[0]) == 8:
      try:
        open(self.files[0],'r').close()
      except:
        print("downloading %s from http://rcsb.org" % self.files[0])
        fetchstr = "http://www.rcsb.org/pdb/cgi/export.cgi/%s?format=PDB&pdbId=2cpk&compression=None" % self.files[0]
        urllib.urlretrieve(fetchstr,self.files[0])

    if self.data and len(self.files): self.read_template(self.files[0])

  # --------------------------------------------------------------------
  # write a single large PDB file for concatenating all input data or files
  # if data exists:
  #   only selected atoms returned by extract
  #   atoms written in order they appear in snapshot
  #   atom only written if its tag is in PDB template file
  # if no data:
  #   concatenate all input files to one output file

  def one(self,*args):
    if len(args) == 0: file = "tmp.pdb"
    elif args[0][-4:] == ".pdb": file = args[0]
    else: file = args[0] + ".pdb"

    f = open(file,'w')

    # use template PDB file with each snapshot

    if self.data:
      n = flag = 0
      while 1:
        which,time,flag = self.data.iterator(flag)
        if flag == -1: break
        self.convert(f,which)
        print("END",file=f)
        print(time,end='')
        sys.stdout.flush()
        n += 1

    else:
      for file in self.files:
        f.write(open(file,'r').read())
        print("END",file=f)
        print(file,end='')
        sys.stdout.flush()

    f.close()
    print("\nwrote %d datasets to %s in PDB format" % (n,file))

  # --------------------------------------------------------------------
  # write series of numbered PDB files
  # if data exists:
  #   only selected atoms returned by extract
  #   atoms written in order they appear in snapshot
  #   atom only written if its tag is in PDB template file
  # if no data:
  #   just copy all input files to output files

  def many(self,*args):
    if len(args) == 0: root = "tmp"
    else: root = args[0]

    if self.data:
      n = flag = 0
      while 1:
        which,time,flag = self.data.iterator(flag)
        if flag == -1: break

        if n < 10:
          file = root + "000" + str(n)
        elif n < 100:
          file = root + "00" + str(n)
        elif n < 1000:
          file = root + "0" + str(n)
        else:
          file = root + str(n)
        file += ".pdb"

        f = open(file,'w')
        self.convert(f,which)
        f.close()

        print(time,end='')
        sys.stdout.flush()
        n += 1

    else:
      n = 0
      for infile in self.files:
        if n < 10:
          file = root + "000" + str(n)
        elif n < 100:
          file = root + "00" + str(n)
        elif n < 1000:
          file = root + "0" + str(n)
        else:
          file = root + str(n)
        file += ".pdb"

        f = open(file,'w')
        f.write(open(infile,'r').read())
        f.close()
        print(file,end='')
        sys.stdout.flush()

        n += 1

    print("\nwrote %d datasets to %s*.pdb in PDB format" % (n,root))

  # --------------------------------------------------------------------
  # write a single PDB file
  # if data exists:
  #   time is timestamp in snapshot
  #   only selected atoms returned by extract
  #   atoms written in order they appear in snapshot
  #   atom only written if its tag is in PDB template file
  # if no data:
  #   time is index into list of input PDB files
  #   just copy one input file to output file

  def single(self,time,*args):
    if len(args) == 0: file = "tmp.pdb"
    elif args[0][-4:] == ".pdb": file = args[0]
    else: file = args[0] + ".pdb"
    f = open(file,'w')

    if self.data:
      which = self.data.findtime(time)
      self.convert(f,which)
    else:
      f.write(open(self.files[time],'r').read())

    f.close()

  # --------------------------------------------------------------------
  # iterate over list of input files or selected snapshots
  # latter is done via data objects iterator

  def iterator(self,flag):
    if not self.data:
      if not flag: self.iterate = 0
      else:
        self.iterate += 1
        if self.iterate > len(self.files): return 0,0,-1
      return self.iterate,self.iterate,1

    return self.data.iterator(flag)

  # --------------------------------------------------------------------
  # read a PDB file and store ATOM lines

  def read_template(self,file):
    lines = open(file,'r').readlines()
    self.atomlines = {}
    for line in lines:
      if line.find("ATOM") == 0:
        tag = int(line[4:11])
        begin = line[:30]
        end = line[54:]
        self.atomlines[tag] = (begin,end)

  # --------------------------------------------------------------------
  # convert one set of atoms to PDB format and write to f

  def convert(self,f,which):
    time,box,atoms,bonds,tris,lines = self.data.viz(which)
    if len(self.files):
      for atom in atoms:
        id = atom[0]
        if id in self.atomlines:
          (begin,end) = self.atomlines[id]
          line = "%s%8.3f%8.3f%8.3f%s" % (begin,atom[2],atom[3],atom[4],end)
          print(line,file=f,end='')
    else:
      for atom in atoms:
        begin = "ATOM %6d %2d   R00     1    " % (atom[0],atom[1])
        middle = "%8.3f%8.3f%8.3f" %  (atom[2],atom[3],atom[4])
        end = "  1.00  0.00    NONE"
        print(begin+middle+end,file=f)
