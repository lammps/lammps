#!/bin/env python

# convert a Tinker XYZ file + PRM file to a LAMMPS data file

# Syntax: python tinker2lmp.py -switch args -switch args ...
#   -xyz file = Tinker XYZ file name (required)
#   -amoeba file = AMOEBA PRM force field file name (required, or hippo)
#   -hippo file = HIPPO PRM force field file name (required, or amoeba)
#   -data file = LAMMPS data file to output (required)
#   -nopbc = non-periodic system (default)
#   -pbc xhi yhi zhi = periodic system from 0 to hi in each dimension (optional)

# Author: Steve Plimpton

import sys,os,math
path = os.environ["LAMMPS_PYTHON_TOOLS"]
sys.path.append(path)
from data import data

BIG = 1.0e20
DELTA = 0.001     # delta on LAMMPS shrink-wrap box size, in Angstroms

# ----------------------
# methods and classes
# ----------------------

# print an error message and exit

def error(txt=""):
  if not txt:
    print "Syntax: tinker2lmp.py -switch args ..."
    print "  -xyz file"
    print "  -amoeba file"
    print "  -hippo file"
    print "  -data file"
    print "  -nopbc"
    print "  -pbc xhi yhi zhi"
  else: print "ERROR:",txt
  #sys.exit()

# read and store values from a Tinker xyz file

class XYZfile:
  def __init__(self,file):
    lines = open(file,'r').readlines()
    natoms = int(lines[0].split()[0])
    id = []
    type = []
    x = []
    y = []
    z = []
    bonds = []
    
    for line in lines[1:natoms+1]:
      words = line.split()
      id.append(int(words[0]))
      x.append(words[2])
      y.append(words[3])
      z.append(words[4])
      type.append(int(words[5]))
      blist = words[6:]
      blist = [int(one) for one in blist]
      bonds.append(blist)
      
    self.natoms = natoms
    self.id = id
    self.type = type
    self.x = x
    self.y = y
    self.z = z
    self.bonds = bonds

  # triplet of atoms in an angle = atom 1,2,3
  # atom2 is center atom
  # nbonds = number of atoms which atom2 is bonded to
  # hcount = # of H atoms which atom2 is bonded to, excluding atom1 and atom3
    
  def angle_hbond_count(self,atom1,atom2,atom3,lmptype,lmpmass):
    bondlist = self.bonds[atom2-1]
    nbonds = len(bondlist)
    
    hcount = 0
    for bondID in bondlist:
      if atom1 == bondID: continue
      if atom3 == bondID: continue
      massone = lmpmass[lmptype[bondID-1]-1]
      if int(massone+0.5) == 1: hcount += 1
      
    return nbonds,hcount
      
# read and store select values from a Tinker force field PRM file
# ntypes = # of Tinker types
# per-type values: class and mass
# scalar force field params in Force Field Definition section
# bond, angle, dihedral coeffs indexed by Tinker classes

class PRMfile:
  def __init__(self,file):
    lines = open(file,'r').readlines()
    self.nlines = len(lines)
    self.lines = lines

    self.force_field_definition()
    self.classes,self.masses = self.peratom()
    self.polgroup = self.polarize()
    self.bondparams = self.bonds()
    self.angleparams = self.angles()
    self.bondangleparams = self.bondangles()
    self.torsionparams = self.torsions()
    self.ureyparams = self.ureybonds()
    self.opbendparams = self.opbend()
    self.ntypes = len(self.masses)

  def force_field_definition(self):
    iline = self.find_section("Force Field Definition")
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "bond-cubic": self.bond_cubic = float(words[1])
        elif words[0] == "bond-quartic": self.bond_quartic = float(words[1])
        elif words[0] == "angle-cubic": self.angle_cubic = float(words[1])
        elif words[0] == "angle-quartic": self.angle_quartic = float(words[1])
        elif words[0] == "angle-pentic": self.angle_pentic = float(words[1])
        elif words[0] == "angle-sextic": self.angle_sextic = float(words[1])
      iline += 1
      
  def peratom(self):
    classes = []
    masses = []
    iline = self.find_section("Atom Type Definitions")
    if iline < 0: return classes,masses
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        # NOTE: assumes atom entries are numbered consecutively
        if words[0] == "atom":
          classes.append(int(words[2]))
          masses.append(float(words[-2]))
      iline += 1
    return classes,masses

  def polarize(self):
    polgroup = []
    iline = self.find_section("Dipole Polarizability Parameters")
    if iline < 0: return polgroup
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        # trim off any end-of-line comments
        if "!!" in words: words = words[:words.index("!!")]
        # NOTE: assumes polarize entries are numbered consecutively
        if words[0] == "polarize":
          if amoeba: bondtypes = words[4:]
          if hippo: bondtypes = words[3:]
          bondtypes = [int(one) for one in bondtypes]
          polgroup.append(bondtypes)
      iline += 1
    return polgroup

  # convert PRMfile params to LAMMPS bond_style class2 params
  
  def bonds(self):
    params = []
    iline = self.find_section("Bond Stretching Parameters")
    if iline < 0: return params
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "bond":
          class1 = int(words[1])
          class2 = int(words[2])
          value1 = float(words[3])
          value2 = float(words[4])
          lmp1 = value2
          lmp2 = value1
          lmp3 = self.bond_cubic * value1
          lmp4 = self.bond_quartic * value1
          params.append((class1,class2,lmp1,lmp2,lmp3,lmp4))
      iline += 1
    return params
    
  # convert PRMfile params to LAMMPS angle_style class2/p6 params
  # line may have prefactor plus 1,2,3 angle0 params
  # save prefactor/angle0 pairs as option 1,2,3
  
  def angles(self):
    r2d = 180.0 / math.pi
    params = []
    iline = self.find_section("Angle Bending Parameters")
    if iline < 0: return params
    iline += 3
    while iline < self.nlines:
      line = self.lines[iline]
      line = line[:line.find("!!")]
      words = line.split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "angle" or words[0] == "anglep":
          if words[0] == "angle": pflag = 0
          if words[0] == "anglep": pflag = 1
          
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          value1 = float(words[4])
          value2 = float(words[5])
          option1 = ()
          option2 = ()
          option3 = ()

          lmp1 = value2
          lmp2 = value1
          lmp3 = self.angle_cubic * value1 * r2d
          lmp4 = self.angle_quartic * value1 * r2d*r2d
          lmp5 = self.angle_pentic * value1 * r2d*r2d*r2d
          lmp6 = self.angle_sextic * value1 * r2d*r2d*r2d*r2d
          
          option1 = (pflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)
          
          if len(words) >= 7:
            value3 = float(words[6])
            lmp1 = value3
            option2 = (pflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)
            
          if len(words) == 8:
            value4 = float(words[7])
            lmp1 = value4
            option3 = (pflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)

          if not option2 and not option3:
            params.append((class1,class2,class3,[option1]))
          elif not option3:
            params.append((class1,class2,class3,[option1,option2]))
          else:
            params.append((class1,class2,class3,[option1,option2,option3]))
          
      iline += 1
    return params

  # convert PRMfile params to LAMMPS angle_style class2/p6 bondangle params
  # lmp3,lmp4 = equilibrium bond lengths for 2 bonds in angle
  #   need to find these values in self.bondparams
  #   put them in a dictionary for efficient searching
  
  def bondangles(self):
    params = []
    iline = self.find_section("Stretch-Bend Parameters")
    if iline < 0: return params
    iline += 3

    bdict = {}
    for m,bparams in enumerate(self.bondparams):
      bdict[(bparams[0],bparams[1])] = bparams[2]
    
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "strbnd":
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          value1 = float(words[4])
          value2 = float(words[5])
          lmp1 = value1
          lmp2 = value2
          
          if (class1,class2) in bdict:
            lmp3 = bdict[(class1,class2)]
          elif (class2,class1) in bdict:
            lmp3 = bdict[(class2,class1)]
          else:
            error("1st bond in BondAngle term not found: %d %d %d" % \
                  (class1,class2,class3))
            # NOTE: just for debugging
            lmp3 = 0.0

          if (class2,class3) in bdict:
            lmp4 = bdict[(class2,class3)]
          elif (class3,class2) in bdict:
            lmp4 = bdict[(class3,class2)]
          else:
            error("2nd bond in BondAngle term not found: %d %d %d" % \
                  (class1,class2,class3))
            # NOTE: just for debugging
            lmp4 = 0.0
            
          params.append((class1,class2,class3,lmp1,lmp2,lmp3,lmp4))
      iline += 1
    return params

  # convert PRMfile params to LAMMPS bond_style class2 bond params
  # coeffs for K2,K3 = 0.0 since Urey-Bradley is simple harmonic
  
  def ureybonds(self):
    params = []
    iline = self.find_section("Urey-Bradley Parameters")
    if iline < 0: return params
    iline += 3

    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "ureybrad":
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          value1 = float(words[4])
          value2 = float(words[5])
          lmp1 = value1
          lmp2 = value2
            
          params.append((class1,class2,class3,lmp1,lmp2,0.0,0.0))
      iline += 1
    return params

  # Dihedral interactions
  
  def torsions(self):
    params = []
    iline = self.find_section("Torsional Parameters")
    if iline < 0: return params
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "torsion":
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          class4 = int(words[4])
          
          if len(words) <= 5:
            error("Torsion has no params: %d %d %d %d" % \
                  (class1,class2,class3,class4))
          if (len(words)-5) % 3: 
            error("Torsion does not have triplets of params: %d %d %d %d" % \
                  (class1,class2,class3,class4))

          mfourier = (len(words)-5) / 3
          oneparams = [class1,class2,class3,class4,mfourier]

          for iset in range(mfourier):
            value1 = float(words[5 + iset*3 + 0])
            value2 = float(words[5 + iset*3 + 1])
            value3 = int(words[5 + iset*3 + 2])
            lmp1 = value1/2.0
            lmp2 = value3
            lmp3 = value2
            oneparams += [lmp1,lmp2,lmp3]
            
          params.append(oneparams)
      iline += 1
    return params

  # Improper or out-of-plane bend interactions
  
  def opbend(self):
    params = []
    iline = self.find_section("Out-of-Plane Bend Parameters")
    if iline < 0: return params
    iline += 3
    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "opbend":
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          class4 = int(words[4])
          value1 = float(words[5])
          lmp1 = value1
          params.append((class1,class2,class3,class4,lmp1))
      iline += 1
    return params
    
  def find_section(self,txt):
    txt = "##  %s  ##" % txt
    for iline,line in enumerate(self.lines):
      if txt in line: return iline
    return -1

# ----------------------------------------
# main program
# ----------------------------------------
    
args = sys.argv[1:]
narg = len(args)

# process args

amoeba = hippo = 0
xyzfile = ""
prmfile = ""
datafile = ""
pbcflag = 0

iarg = 0
while iarg < narg:
  if args[iarg] == "-xyz":
    if iarg + 2 > narg: error()
    xyzfile = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-amoeba":
    if iarg + 2 > narg: error()
    amoeba = 1
    prmfile = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-hippo":
    if iarg + 2 > narg: error()
    hippo = 1
    prmfile = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-data":
    if iarg + 2 > narg: error()
    datafile = args[iarg+1]
    iarg += 2
  elif args[iarg] == "-nopbc":
    pbcflag = 0
    iarg += 1
  elif args[iarg] == "-pbc":
    if iarg + 4 > narg: error()
    pbcflag = 1
    xhi = float(args[iarg+1])
    yhi = float(args[iarg+2])
    zhi = float(args[iarg+3])
    boxhi = [xhi,yhi,zhi]
    iarg += 4
  else: error()

# error check

if not xyzfile: error("xyzfile not specified")
if not prmfile: error("prmfile not specified")
if not datafile: error("datafile not specified")

# read Tinker xyz and prm files

xyz = XYZfile(xyzfile)
prm = PRMfile(prmfile)

# create LAMMPS box bounds based on pbcflag

natoms = xyz.natoms
x = xyz.x
y = xyz.y
z = xyz.z

if pbcflag:
  boxlo = [0,0,0]
else:
  xlo = ylo = zlo = BIG
  xhi = yhi = zhi = -BIG
  for i in xrange(natoms):
    xlo = min(xlo,float(x[i]))
    ylo = min(ylo,float(y[i]))
    zlo = min(zlo,float(z[i]))
    xhi = max(xhi,float(x[i]))
    yhi = max(yhi,float(y[i]))
    zhi = max(zhi,float(z[i]))  
  boxlo = [xlo-DELTA,ylo-DELTA,zlo-DELTA]
  boxhi = [xhi+DELTA,yhi+DELTA,zhi+DELTA]

# ----------------------------------------
# create LAMMPS atom types for each unique Tinker per-type mass
# NOTE: maybe should assign LAMMPS types in a different way,
#       e.g. one LAMMPS type for each used Tinker type
# ----------------------------------------

# ntypes = # of LAMMPS atoms types = unique Tinker masses
# lmptype = which LAMMPS type for each atom (1 to Ntypes)
# lmpmass = list of per-type masses
# ttype = list of Tinker types for each atom (1 to prm.ntypes)
# tink2lmp = mapping of Tinker types to LAMMPS types

natoms = xyz.natoms
ttype = xyz.type
tmasses = prm.masses

ntypes = 0
lmptype = []
lmpmass = []
tink2lmp = {}

for itype in ttype:
  if itype not in tink2lmp:
    mass = tmasses[itype-1]
    if mass not in lmpmass:
      ntypes += 1
      lmpmass.append(mass)
      jtype = ntypes
    else: jtype = lmpmass.index(mass) + 1
    tink2lmp[itype] = jtype
  lmptype.append(tink2lmp[itype])

# ----------------------------------------
# identify molecules from Tinker bond connectivity
# ----------------------------------------

# molID = which molecule each atom is in (1 to Nmol)
# use stack to store IDs of atoms to recursively loop over

natoms = xyz.natoms
id = xyz.id
bonds = xyz.bonds

molID = natoms*[0]
nmol = 0

for i in id:
  if molID[i-1] > 0: continue
  nmol += 1
  molID[i-1] = nmol
  stack = [i]
  while stack:
    j = stack.pop()
    for k in bonds[j-1]:
      if molID[k-1] == 0:
        molID[k-1] = nmol
        stack.append(k)

# ----------------------------------------
# create lists of bonds, angles, dihedrals, impropers
# ----------------------------------------

# create blist = list of bonds
# avoid double counting by requiring atom1 < atom2

id = xyz.id
type = xyz.type
bonds = xyz.bonds

blist = []

for atom1 in id:
  for atom2 in bonds[atom1-1]:
    if atom1 < atom2:
      blist.append((atom1,atom2))

# create alist = list of angles
# generate topology by double loop over bonds of center atom2
# avoid double counting by requiring atom1 < atom3

id = xyz.id
type = xyz.type
bonds = xyz.bonds

alist = []

for atom2 in id:
  for atom1 in bonds[atom2-1]:
    for atom3 in bonds[atom2-1]:
      if atom3 == atom1: continue
      if atom1 < atom3:
        alist.append((atom1,atom2,atom3))

# create dlist = list of dihedrals
# generate topology via triple loop over neighbors of dihedral atom2
#   double loop over bonds of atom2
#   additional loop over bonds of atom3
# avoid double counting the reverse dihedral by use of ddict dictionary

id = xyz.id
type = xyz.type
bonds = xyz.bonds

dlist = []
ddict = {}

for atom2 in id:
  for atom1 in bonds[atom2-1]:
    for atom3 in bonds[atom2-1]:
      if atom3 == atom1: continue
      for atom4 in bonds[atom3-1]:
        if atom4 == atom2 or atom4 == atom1: continue
        if (atom4,atom3,atom2,atom1) in ddict: continue
        dlist.append((atom1,atom2,atom3,atom4))
        ddict[(atom1,atom2,atom3,atom4)] = 1

# create olist = list of out-of-plane impropers
# generate topology by triple loop over bonds of center atom2
# atom2 must have 3 or more bonds to be part of an improper
# avoid double counting by requiring atom3 < atom4
# this is since in Tinker the final 2 atoms in the improper are interchangeable

id = xyz.id
type = xyz.type
bonds = xyz.bonds

olist = []

for atom2 in id:
  if len(bonds[atom2-1]) < 3: continue
  for atom1 in bonds[atom2-1]:
    for atom3 in bonds[atom2-1]:
      for atom4 in bonds[atom2-1]:
        if atom1 == atom3: continue
        if atom1 == atom4: continue
        if atom3 >= atom4: continue
        olist.append((atom1,atom2,atom3,atom4))

# ----------------------------------------
# create list of Urey-Bradley triplet matches
# ----------------------------------------

# scan list of angles to find triplets that match UB parameters
# if match, add it to UB bond list

type = xyz.type
classes = prm.classes

ublist = []

ubdict = {}
for m,params in enumerate(prm.ureyparams):
  ubdict[(params[0],params[1],params[2])] = (m,params)

for atom1,atom2,atom3 in alist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]

  if (c1,c2,c3) in ubdict:
    ublist.append((atom1,atom2,atom3))
  elif (c3,c2,c1) in ubdict:
    ublist.append((atom3,atom2,atom1))

# ----------------------------------------
# create lists of bond/angle/dihedral/improper types
# ----------------------------------------

# generate btype = LAMMPS type of each bond
# generate bparams = LAMMPS params for each bond type
# flags[i] = which LAMMPS bond type (1-N) the Ith Tinker PRM file bond is
#        0 = none
# convert prm.bondparams to a dictionary for efficient searching
# key = (class1,class2)
# value = (M,params) where M is index into prm.bondparams

id = xyz.id
type = xyz.type
classes = prm.classes

bdict = {}
for m,params in enumerate(prm.bondparams):
  bdict[(params[0],params[1])] = (m,params)
         
flags = len(prm.bondparams)*[0]
btype = []
bparams = []

for atom1,atom2 in blist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  
  if (c1,c2) in bdict: m,params = bdict[(c1,c2)]
  elif (c2,c1) in bdict: m,params = bdict[(c2,c1)]
  else: error("Bond not found: %d %d: %d %d" % (atom1,atom2,c1,c2))
    
  if not flags[m]:
    v1,v2,v3,v4 = params[2:]
    bparams.append((v1,v2,v3,v4))
    flags[m] = len(bparams)
  btype.append(flags[m])

# extend btype and bparams with Urey-Bradley H-H bond info
# loop over ublist of UB bonds
# convert prm.ureyparams to a dictionary for efficient searching
# key = (class1,class2,class3)
# value = (M,params) where M is index into prm.ureyparams
# also add the c1,c3 H-H bond to blist

id = xyz.id
type = xyz.type
classes = prm.classes

ubdict = {}
for m,params in enumerate(prm.ureyparams):
  ubdict[(params[0],params[1],params[2])] = (m,params)

flags = len(prm.ureyparams)*[0]

for atom1,atom2,atom3 in ublist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]
  
  if (c1,c2,c3) in ubdict:
    m,params = ubdict[(c1,c2,c3)]
    blist.append((c1,c3))
  elif (c3,c2,c1) in ubdict:
    m,params = ubdict[(c3,c2,c1)]
    blist.append((c3,c1))
    
  if not flags[m]:
    v1,v2,v3,v4 = params[3:]
    bparams.append((v1,v2,v3,v4))
    flags[m] = len(bparams)
  btype.append(flags[m])
  
# generate atype = LAMMPS type of each angle
# generate aparams = LAMMPS params for each angle type
# flags[i] = which LAMMPS angle type (1-N) the Tinker FF file angle I is
#        0 = none
# Tinker FF file angle entries can have 1, 2, or 3 options
# noptions = total # of Tinker FF file entries with options included
# convert prm.angleparams to a dictionary for efficient searching
# key = (class1,class2)
# value = (M,params) where M is index into prm.angleparams

id = xyz.id
type = xyz.type
classes = prm.classes

adict = {}
noptions = 0
for m,params in enumerate(prm.angleparams):
  adict[(params[0],params[1],params[2])] = (noptions,params)
  n = len(params[3])
  noptions += n
  
flags = noptions*[0]
atype = []
aparams = []

for i,one in enumerate(alist):
  atom1,atom2,atom3 = one
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]

  if (c1,c2,c3) in adict or (c3,c2,c1) in adict:
    if (c1,c2,c3) in adict: m,params = adict[(c1,c2,c3)]

    # IMPORTANT subtlety
    # flip order of 3 atoms in alist if the angle
    #   matches Angle Bending section of PRM file in reverse order
    # necessary b/c BondAngle coeffs will be generated with r1,r2 params
    #   from Bond Stretching section of PRM file
    # since in general r1 != r2, the LAMMPS AngleAmoeba class requires
    #   the 3 atoms in the angle be in the order that matches r1 and r2
    
    if (c3,c2,c1) in adict:
      m,params = adict[(c3,c2,c1)]
      alist[i] = (atom3,atom2,atom1)
      
    # params is a sequence of 1 or 2 or 3 options
    # which = which of 1,2,3 options this atom triplet matches
    # for which = 2 or 3, increment m to index correct position in flags
    # how match is determined:
    #   if 2 options:
    #     require atom2 have 3 bond partners, including atom1 and atom3
    #     option 1 if additional bond is not to an H atom
    #     option 2 if additional bond is to an H atom
    #   if 3 options:
    #     require atom2 have 4 bond partners, including atom1 and atom3
    #     option 1 if neither of 2 additional bonds is to an H atom
    #     option 2 if one of 2 additional bonds is to an H atom
    #     option 3 if both of 2 additional bonds is to an H atom

    if len(params[3]) == 1:
      which = 1
      
    elif len(params[3]) == 2:
      nbonds,hcount = xyz.angle_hbond_count(atom1,atom2,atom3,lmptype,lmpmass)
      
      if nbonds != 3: 
        print "Center angle atom has wrong bond count"
        print "  angle atom IDs:",atom1,atom2,atom3
        print "  angle atom classes:",c1,c2,c3
        print "  Tinker FF file param options:",len(params[3])
        print "  Nbonds and hydrogen count:",nbonds,hcount
        #sys.exit()      NOTE: allow this for now

      if hcount == 0: which = 1
      elif hcount == 1:
        which = 2
        m += 1

      print "3-bond angle"
      print "  angle atom IDs:",atom1,atom2,atom3
      print "  angle atom classes:",c1,c2,c3
      print "  Tinker FF file param options:",len(params[3])
      print "  Nbonds and hydrogen count:",nbonds,hcount
      print "  which:",which,m

    elif len(params[3]) == 3:
      nbonds,hcount = xyz.angle_hbond_count(atom1,atom2,atom3,lmptype,lmpmass)
      
      if nbonds != 4: 
        print "Center angle atom has wrong bond count"
        print "  angle atom IDs:",atom1,atom2,atom3
        print "  angle atom classes:",c1,c2,c3
        print "  Tinker FF file param options:",len(params[3])
        print "  Nbonds and hydrogen count:",nbonds,hcount
        #sys.exit()     NOTE: allow this for now
        
      if hcount == 0: which = 1
      elif hcount == 1:
        which = 2
        m += 1
      elif hcount == 2:
        which = 3
        m += 2

  else:
    error("Angle not found: %d %d %d: %d %d %d" % (atom1,atom2,atom3,c1,c2,c3))
    
  if not flags[m]:
    pflag,v1,v2,v3,v4,v5,v6 = params[3][which-1]
    aparams.append((pflag,v1,v2,v3,v4,v5,v6))
    flags[m] = len(aparams)
  atype.append(flags[m])

# augment the aparams with bond-angle cross terms from bondangleparams
# generate baparams = LAMMPS bond-angle params for each angle type
# sbdict = dictionary for angle tuples in bongangleparams

sbdict = {}
for v1,v2,v3,v4,v5,v6,v7 in prm.bondangleparams:
  if (v1,v2,v3) in sbdict: continue
  sbdict[(v1,v2,v3)] = (v4,v5,v6,v7)

baparams = []

for itype in range(len(aparams)):
  iangle = atype.index(itype+1) 
  atom1,atom2,atom3 = alist[iangle]
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]

  if (c1,c2,c3) in sbdict:
    n1,n2,r1,r2 = sbdict[(c1,c2,c3)]
  elif (c3,c2,c1) in sbdict:
    n1,n2,r1,r2 = sbdict[(c3,c2,c1)]
  else:
    print "Bond-stretch angle triplet not found: %d %d %d" % (c1,c2,c3)
    n1,n2,r1,r2 = 4*[0.0]
    
  baparams.append((n1,n2,r1,r2))

# generate dtype = LAMMPS type of each dihedral
# generate dparams = LAMMPS params for each dihedral type
# flags[i] = which LAMMPS dihedral type (1-N) the Tinker FF file dihedral I is
#        0 = none
# convert prm.torsionparams to a dictionary for efficient searching
# key = (class1,class2,class3,class4)
# value = (M,params) where M is index into prm.torsionparams

id = xyz.id
type = xyz.type
classes = prm.classes

ddict = {}
for m,params in enumerate(prm.torsionparams):
  ddict[(params[0],params[1],params[2],params[3])] = (m,params)

flags = len(prm.torsionparams)*[0]
dtype = []
dparams = []

for atom1,atom2,atom3,atom4 in dlist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  type4 = type[atom4-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]
  c4 = classes[type4-1]
  
  if (c1,c2,c3,c4) in ddict: m,params = ddict[(c1,c2,c3,c4)]
  elif (c4,c3,c2,c1) in ddict: m,params = ddict[(c4,c3,c2,c1)]
  else:
    error("Dihedral not found: %d %d %d %d: %d %d %d %d" % \
          (atom1,atom2,atom3,atom4,c1,c2,c3,c4))

  if not flags[m]:
    oneparams = params[4:]
    dparams.append(oneparams)
    flags[m] = len(dparams)
  dtype.append(flags[m])

# generate otype = LAMMPS type of each out-of-plane improper
# generate oparams = LAMMPS params for each improper type
# flags[i] = which LAMMPS improper type (1-N) the Tinker FF file improper I is
#        0 = none
# convert prm.opbendparams to a dictionary for efficient searching
# key = (class1,class2)
# value = (M,params) where M is index into prm.opbendparams

id = xyz.id
type = xyz.type
classes = prm.classes

odict = {}
for m,params in enumerate(prm.opbendparams):
  odict[(params[0],params[1])] = (m,params)

flags = len(prm.opbendparams)*[0]
otype = []
oparams = []
olist_reduced = []

for atom1,atom2,atom3,atom4 in olist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  type4 = type[atom4-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]
  c4 = classes[type4-1]

  # 4-tuple is only an improper if matches an entry in PRM file
  # olist_reduced = list of just these 4-tuples
  
  if (c1,c2) in odict:
    m,params = odict[(c1,c2)]
    olist_reduced.append((atom1,atom2,atom3,atom4))
    
    if not flags[m]:
      oneparams = params[4:]
      oparams.append(oneparams)
      flags[m] = len(oparams)
    otype.append(flags[m])

# replace original olist with reduced version

olist = olist_reduced
    
# ----------------------------------------
# assign each atom to a Tinker group
# NOTE: doing this inside LAMMPS now
#       comment out unless need to test
# ----------------------------------------

# ngroup = # of groups
# tgroup[i] = groupID (1 to Ngroups) for each atom
# use stack to store IDs of atoms to recursively loop over
# do not set tgroup for an atom if already assigned to a group
# only add atoms J (bonded to M) to stack
#   whose Tinker type is in polgroup list for atom M's Tinker type

#natoms = xyz.natoms
#id = xyz.id
#bonds = xyz.bonds
#ttype = xyz.type
#polgroup = prm.polgroup

#tgroup = natoms*[0]
#ngroups = 0

#for i in id:
#  if tgroup[i-1] > 0: continue
#
#  ngroups += 1
#  groupID = ngroups
#  stack = [i]
#  while stack:
#    m = stack.pop()
#    if tgroup[m-1] > 0: continue
#    tgroup[m-1] = groupID
#    for j in bonds[m-1]:
#      if tgroup[j-1] > 0: continue
#      if ttype[j-1] not in polgroup[ttype[m-1]-1]: continue
#      stack.append(j)
      
# ----------------------------------------
# write LAMMPS data file via Pizza.py data class
# ----------------------------------------

d = data()

natoms = xyz.natoms
id = xyz.id
x = xyz.x
y = xyz.y
z = xyz.z
ttype = xyz.type

nbonds = len(blist)
nangles = len(alist)
ndihedrals = len(dlist)
nimpropers = len(olist)

# data file header values

d.title = "LAMMPS data file created from Tinker %s and %s files\n" % \
          (xyzfile,prmfile)

d.headers["atoms"] = natoms
d.headers["atom types"] = ntypes

d.headers["xlo xhi"] = (boxlo[0],boxhi[0])
d.headers["ylo yhi"] = (boxlo[1],boxhi[1])
d.headers["zlo zhi"] = (boxlo[2],boxhi[2])

# data file sections

lines = []
for i,mass in enumerate(lmpmass):
  line = "%d %s" % (i+1,mass)
  lines.append(line+'\n')
d.sections["Masses"] = lines

lines = []
for i,one in enumerate(id):
  line = "%d %d %d %g %s %s %s" % (one,molID[i],lmptype[i],0.0,x[i],y[i],z[i])
  lines.append(line+'\n')
d.sections["Atoms"] = lines

lines = []
for i,one in enumerate(ttype):
  # comment out inclusion of Tinker group, now done by LAMMPS
  #line = "%d %d %d" % (id[i],one,tgroup[i])
  line = "%d %d" % (id[i],one)
  lines.append(line+'\n')
d.sections["Tinker Types"] = lines

if nbonds:
  d.headers["bonds"] = len(blist)
  d.headers["bond types"] = len(bparams)
  
  lines = []
  for i,one in enumerate(bparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["Bond Coeffs"] = lines

  lines = [] 
  for i,one in enumerate(blist):
    line = "%d %d %d %d" % (i+1,btype[i],one[0],one[1])
    lines.append(line+'\n')
  d.sections["Bonds"] = lines

if nangles:
  d.headers["angles"] = len(alist)
  d.headers["angle types"] = len(aparams)
  
  lines = []
  for i,one in enumerate(aparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["Angle Coeffs"] = lines

  lines = []
  for i,one in enumerate(baparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["BondAngle Coeffs"] = lines

  lines = [] 
  for i,one in enumerate(alist):
    line = "%d %d %d %d %d" % (i+1,atype[i],one[0],one[1],one[2])
    lines.append(line+'\n')
  d.sections["Angles"] = lines

if ndihedrals:
  d.headers["dihedrals"] = len(dlist)
  d.headers["dihedral types"] = len(dparams)
  
  lines = []
  for i,one in enumerate(dparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["Dihedral Coeffs"] = lines

  lines = [] 
  for i,one in enumerate(dlist):
    line = "%d %d %d %d %d %d" % (i+1,dtype[i],one[0],one[1],one[2],one[3])
    lines.append(line+'\n')
  d.sections["Dihedrals"] = lines

if nimpropers:
  d.headers["impropers"] = len(olist)
  d.headers["improper types"] = len(oparams)
  
  lines = []
  for i,one in enumerate(oparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["Improper Coeffs"] = lines

  lines = [] 
  for i,one in enumerate(olist):
    line = "%d %d %d %d %d %d" % (i+1,otype[i],one[0],one[1],one[2],one[3])
    lines.append(line+'\n')
  d.sections["Impropers"] = lines

d.write(datafile)

# print stats to screen

print "Natoms =",natoms
print "Ntypes =",ntypes
print "Tinker XYZ types =",len(tink2lmp)
print "Tinker PRM types =",prm.ntypes
#print "Tinker groups =",ngroups
print "Nmol =",nmol
print "Nbonds =",len(blist)
print "Nangles =",len(alist)
print "Ndihedrals =",len(dlist)
print "Nimpropers =",len(olist)
print "Nbondtypes =",len(bparams)
print "Nangletypes =",len(aparams)
print "Ndihedraltypes =",len(dparams)
print "Nimpropertypes =",len(oparams)
