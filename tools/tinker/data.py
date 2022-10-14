# Pizza.py toolkit, https://lammps.github.io/pizza
# Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
#
# Copyright (2005) Sandia Corporation.  Under the terms of Contract
# DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
# certain rights in this software.  This software is distributed under
# the GNU General Public License.

# data tool

from __future__ import print_function
oneline = "Read, write, manipulate LAMMPS data files"

docstr = """
d = data("data.poly")            read a LAMMPS data file, can be gzipped
d = data()           create an empty data file

d.map(1,"id",3,"x")              assign names to atom columns (1-N)

coeffs = d.get("Pair Coeffs")    extract info from data file section
q = d.get("Atoms",4)

  1 arg = all columns returned as 2d array of floats
  2 args = Nth column returned as vector of floats

d.reorder("Atoms",1,3,2,4,5)     reorder columns (1-N) in a data file section

  1,3,2,4,5 = new order of previous columns, can delete columns this way

d.title = "My LAMMPS data file"  set title of the data file
d.headers["atoms"] = 1500        set a header value
d.sections["Bonds"] = lines      set a section to list of lines (with newlines)
d.delete("bonds")        delete a keyword or section of data file
d.delete("Bonds")
d.replace("Atoms",5,vec)         replace Nth column of section with vector
d.newxyz(dmp,1000)       replace xyz in Atoms with xyz of snapshot N

  newxyz assumes id,x,y,z are defined in both data and dump files
    also replaces ix,iy,iz if they are defined

index,time,flag = d.iterator(0/1)          loop over single data file snapshot
time,box,atoms,bonds,tris,lines = d.viz(index)   return list of viz objects

  iterator() and viz() are compatible with equivalent dump calls
  iterator() called with arg = 0 first time, with arg = 1 on subsequent calls
    index = timestep index within dump object (only 0 for data file)
    time = timestep value (only 0 for data file)
    flag = -1 when iteration is done, 1 otherwise
  viz() returns info for specified timestep index (must be 0)
    time = 0
    box = [xlo,ylo,zlo,xhi,yhi,zhi]
    atoms = id,type,x,y,z for each atom as 2d array
    bonds = id,type,x1,y1,z1,x2,y2,z2,t1,t2 for each bond as 2d array
      NULL if bonds do not exist
    tris = NULL
    lines = NULL

d.write("data.new")             write a LAMMPS data file
"""

# History
#   8/05, Steve Plimpton (SNL): original version
#   11/07, added triclinic box support

# ToDo list

# Variables
#   title = 1st line of data file
#   names = dictionary with atom attributes as keys, col #s as values
#   headers = dictionary with header name as key, value or tuple as values
#   sections = dictionary with section name as key, array of lines as values
#   nselect = 1 = # of snapshots

# Imports and external programs

from os import popen

try: tmp = PIZZA_GUNZIP
except: PIZZA_GUNZIP = "gunzip"

# Class definition

class data(object):

  # --------------------------------------------------------------------

  def __init__(self,*list):
    self.nselect = 1

    if len(list) == 0:
      self.title = "LAMMPS data file"
      self.names = {}
      self.headers = {}
      self.sections = {}
      return

    file = list[0]
    if file[-3:] == ".gz": f = popen("%s -c %s" % (PIZZA_GUNZIP,file),'r')
    else: f = open(file)

    self.title = f.readline()
    self.names = {}

    headers = {}
    while 1:
      line = f.readline()
      line = line.strip()
      if len(line) == 0:
        continue
      found = 0
      for keyword in hkeywords:
        if line.find(keyword) >= 0:
          found = 1
          words = line.split()
          if keyword == "xlo xhi" or keyword == "ylo yhi" or \
            keyword == "zlo zhi":
            headers[keyword] = (float(words[0]),float(words[1]))
          elif keyword == "xy xz yz":
            headers[keyword] = \
              (float(words[0]),float(words[1]),float(words[2]))
          else:
            headers[keyword] = int(words[0])
      if not found:
        break

    sections = {}
    while 1:
      found = 0
      for pair in skeywords:
        keyword,length = pair[0],pair[1]
        if keyword == line:
          found = 1
          if length not in headers:
            raise (Exception, "data section %s has no matching header value" % line)
          f.readline()
          list = []
          for i in range(headers[length]): list.append(f.readline())
          sections[keyword] = list
      if not found:
        raise (Exception,"invalid section %s in data file" % line)
      f.readline()
      line = f.readline()
      if not line:
        break
      line = line.strip()

    f.close()
    self.headers = headers
    self.sections = sections

  # --------------------------------------------------------------------
  # assign names to atom columns

  def map(self,*pairs):
    if len(pairs) % 2 != 0:
      raise Exception("data map() requires pairs of mappings")
    for i in range(0,len(pairs),2):
      j = i + 1
      self.names[pairs[j]] = pairs[i]-1

  # --------------------------------------------------------------------
  # extract info from data file fields

  def get(self,*list):
    if len(list) == 1:
      field = list[0]
      array = []
      lines = self.sections[field]
      for line in lines:
        words = line.split()
        values = list(map(float,words))
        array.append(values)
      return array
    elif len(list) == 2:
      field = list[0]
      n = list[1] - 1
      vec = []
      lines = self.sections[field]
      for line in lines:
        words = line.split()
        vec.append(float(words[n]))
      return vec
    else:
      raise Exception("invalid arguments for data.get()")

  # --------------------------------------------------------------------
  # reorder columns in a data file field

  def reorder(self,name,*order):
    n = len(order)
    natoms = len(self.sections[name])
    oldlines = self.sections[name]
    newlines = natoms*[""]
    for index in order:
      for i in range(len(newlines)):
        words = oldlines[i].split()
        newlines[i] += words[index-1] + " "
    for i in range(len(newlines)):
      newlines[i] += "\n"
    self.sections[name] = newlines

  # --------------------------------------------------------------------
  # replace a column of named section with vector of values

  def replace(self,name,icol,vector):
    lines = self.sections[name]
    newlines = []
    j = icol - 1
    for i in range(len(lines)):
      line = lines[i]
      words = line.split()
      words[j] = str(vector[i])
      newline = ' '.join(words) + '\n'
      newlines.append(newline)
    self.sections[name] = newlines

  # --------------------------------------------------------------------
  # replace x,y,z in Atoms with x,y,z values from snapshot ntime of dump object
  # assumes id,x,y,z are defined in both data and dump files
  # also replaces ix,iy,iz if they are defined

  def newxyz(self,dm,ntime):
    nsnap = dm.findtime(ntime)

    dm.sort(ntime)
    x,y,z = dm.vecs(ntime,"x","y","z")

    self.replace("Atoms",self.names['x']+1,x)
    self.replace("Atoms",self.names['y']+1,y)
    self.replace("Atoms",self.names['z']+1,z)

    if "ix" in dm.names and "ix" in self.names:
      ix,iy,iz = dm.vecs(ntime,"ix","iy","iz")
      self.replace("Atoms",self.names['ix']+1,ix)
      self.replace("Atoms",self.names['iy']+1,iy)
      self.replace("Atoms",self.names['iz']+1,iz)

  # --------------------------------------------------------------------
  # delete header value or section from data file

  def delete(self,keyword):

    if keyword in self.headers: del self.headers[keyword]
    elif keyword in self.sections: del self.sections[keyword]
    else: raise Exception("keyword not found in data object")

  # --------------------------------------------------------------------
  # write out a LAMMPS data file

  def write(self,file):
    f = open(file,"w")
    print(self.title, file=f)

    # write any keywords in standard list hkeywords
    #   in the order they are in hkeywords
    # then write any extra keywords at end of header section

    for keyword in hkeywords:
      if keyword in self.headers:
        if keyword == "xlo xhi" or keyword == "ylo yhi" or \
               keyword == "zlo zhi":
          pair = self.headers[keyword]
          print(pair[0],pair[1],keyword, file=f)
        elif keyword == "xy xz yz":
          triple = self.headers[keyword]
          print(triple[0],triple[1],triple[2],keyword, file=f)
        else:
          print(self.headers[keyword],keyword, file=f)

    for keyword in list(self.headers.keys()):
      if keyword not in hkeywords:
        print(self.headers[keyword],keyword, file=f)

    # write any sections in standard list skeywords
    #   in the order they are in skeywords
    # then write any extra sections at end of file

    for pair in skeywords:
      keyword = pair[0]
      if keyword in self.sections:
        print("\n%s\n" % keyword, file=f)
        for line in self.sections[keyword]:
          print(line, end='', file=f)

    skeyfirst = [pair[0] for pair in skeywords]

    for keyword in list(self.sections.keys()):
      if keyword not in skeyfirst:
        print("\n%s\n" % keyword, file=f)
        for line in self.sections[keyword]:
          print(line, end='', file=f)

    f.close()

  # --------------------------------------------------------------------
  # iterator called from other tools

  def iterator(self,flag):
    if flag == 0: return 0,0,1
    return 0,0,-1

  # --------------------------------------------------------------------
  # time query from other tools

  def findtime(self,n):
    if n == 0: return 0
    raise(Exception, "no step %d exists" % (n))

  # --------------------------------------------------------------------
  # return list of atoms and bonds to viz for data object

  def viz(self,isnap):
    if isnap: raise Exception("cannot call data.viz() with isnap != 0")

    id = self.names["id"]
    type = self.names["type"]
    x = self.names["x"]
    y = self.names["y"]
    z = self.names["z"]

    xlohi = self.headers["xlo xhi"]
    ylohi = self.headers["ylo yhi"]
    zlohi = self.headers["zlo zhi"]
    box = [xlohi[0],ylohi[0],zlohi[0],xlohi[1],ylohi[1],zlohi[1]]

    # create atom list needed by viz from id,type,x,y,z

    atoms = []
    atomlines = self.sections["Atoms"]
    for line in atomlines:
      words = line.split()
      atoms.append([int(words[id]),int(words[type]),
                    float(words[x]),float(words[y]),float(words[z])])

    # create list of current bond coords from list of bonds
    # assumes atoms are sorted so can lookup up the 2 atoms in each bond

    bonds = []
    if "Bonds" in self.sections:
      bondlines = self.sections["Bonds"]
      for line in bondlines:
        words = line.split()
        bid,btype   = int(words[0]),int(words[1])
        atom1,atom2 = int(words[2]),int(words[3])
        atom1words  = atomlines[atom1-1].split()
        atom2words  = atomlines[atom2-1].split()
        bonds.append([bid,btype,
                      float(atom1words[x]),float(atom1words[y]),
                      float(atom1words[z]),
                      float(atom2words[x]),float(atom2words[y]),
                      float(atom2words[z]),
                      float(atom1words[type]),float(atom2words[type])])

    tris = []
    lines = []
    return 0,box,atoms,bonds,tris,lines

  # --------------------------------------------------------------------
  # return box size

  def maxbox(self):
    xlohi = self.headers["xlo xhi"]
    ylohi = self.headers["ylo yhi"]
    zlohi = self.headers["zlo zhi"]
    return [xlohi[0],ylohi[0],zlohi[0],xlohi[1],ylohi[1],zlohi[1]]

  # --------------------------------------------------------------------
  # return number of atom types

  def maxtype(self):
    return self.headers["atom types"]

# --------------------------------------------------------------------
# standard data file keywords, both header and main sections

hkeywords = ["atoms","ellipsoids","lines","triangles","bodies",
             "bonds","angles","dihedrals","impropers",
             "atom types","bond types","angle types","dihedral types",
             "improper types",
             "xlo xhi","ylo yhi","zlo zhi","xy xz yz"]

skeywords = [["Masses","atom types"],
             ["Atoms","atoms"],["Ellipsoids","ellipsoids"],
             ["Lines","lines"],["Triangles","triangles"],["Bodies","bodies"],
             ["Velocities","atoms"],
             ["Bonds","bonds"],
             ["Angles","angles"],
             ["Dihedrals","dihedrals"],
             ["Impropers","impropers"],
             ["Pair Coeffs","atom types"],
             ["Bond Coeffs","bond types"],
             ["Angle Coeffs","angle types"],
             ["Dihedral Coeffs","dihedral types"],
             ["Improper Coeffs","improper types"],
             ["BondBond Coeffs","angle types"],
             ["BondAngle Coeffs","angle types"],
             ["MiddleBondTorsion Coeffs","dihedral types"],
             ["EndBondTorsion Coeffs","dihedral types"],
             ["AngleTorsion Coeffs","dihedral types"],
             ["AngleAngleTorsion Coeffs","dihedral types"],
             ["BondBond13 Coeffs","dihedral types"],
             ["AngleAngle Coeffs","improper types"]]
