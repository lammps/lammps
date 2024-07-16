#!/bin/env python

# convert a Tinker XYZ file + PRM file to a LAMMPS data file

# Syntax: python tinker2lmp.py -switch args -switch args ...
#   -xyz file = Tinker XYZ file name (required)
#   -amoeba file = AMOEBA PRM force field file name (required, or hippo)
#   -hippo file = HIPPO PRM force field file name (required, or amoeba)
#   -data file = LAMMPS data file to output (required)
#   -bitorsion file = LAMMPS fix bitorsion file to output (required if BiTorsions)
#   -nopbc = non-periodic system (default)
#   -pbc xhi yhi zhi = periodic system from 0 to hi in each dimension (optional)
#   -rep Nx Ny Nz outfile = replicate periodic system by Nx by Ny by Nz (optional)
#      outfile = new xyz file to write out, Null if no file

# Author: Steve Plimpton

from __future__ import print_function
import sys,os,math
from data import data

BIG = 1.0e20
DELTA = 0.001     # delta on LAMMPS shrink-wrap box size, in Angstroms

# ----------------------
# methods and classes
# ----------------------

# print an error message and exit

def error(txt="""
Syntax: tinker2lmp.py -switch args ...
  -xyz file
  -amoeba file
  -hippo file
  -data file
  -bitorsion file
  -nopbc
  -pbc xhi yhi zhi"""):
  sys.exit("ERROR: " + txt)

# read and store values from a Tinker xyz file

class XYZfile(object):
  def __init__(self,file):
    lines = open(file,'r').readlines()
    header = lines[0]
    natoms = int(lines[0].split()[0])
    id = []
    label = []
    type = []
    x = []
    y = []
    z = []
    bonds = []

    for line in lines[1:natoms+1]:
      words = line.split()
      id.append(int(words[0]))
      label.append(words[1])
      x.append(words[2])
      y.append(words[3])
      z.append(words[4])
      type.append(int(words[5]))
      blist = words[6:]
      blist = [int(one) for one in blist]
      bonds.append(blist)

    self.header = header
    self.natoms = natoms
    self.id = id
    self.label = label
    self.type = type
    self.x = x
    self.y = y
    self.z = z
    self.bonds = bonds

  # set bond images flags in each dim for each bond of each atom
  # used for replication of a periodic system
  # 0 = bond partner is closer then half box in dim
  # -1 = bond partner is on other side of lower box bound in dim
  # +1 = bond partner is on other side of upper box bound in dim

  def bond_images(self):

    xhalfsq = 0.25 * boxhi[0]*boxhi[0]
    yhalfsq = 0.25 * boxhi[1]*boxhi[1]
    zhalfsq = 0.25 * boxhi[2]*boxhi[2]

    x = self.x
    y = self.y
    z = self.z
    bonds = self.bonds
    imageflags = []

    for i in range(self.natoms):
      xi = float(x[i])
      yi = float(y[i])
      zi = float(z[i])
      oneflags = []

      for j in bonds[i]:
        xj = float(x[j-1])
        yj = float(y[j-1])
        zj = float(z[j-1])
        dx = xi - xj
        dy = yi - yj
        dz = zi - zj

        if dx*dx <= xhalfsq: ximage = 0
        elif xi < xj: ximage = -1
        elif xi > xj: ximage = 1
        if dy*dy <= yhalfsq: yimage = 0
        elif yi < yj: yimage = -1
        elif yi > yj: yimage = 1
        if dz*dz <= zhalfsq: zimage = 0
        elif zi < zj: zimage = -1
        elif zi > zj: zimage = 1
        oneflags.append((ximage,yimage,zimage))

      imageflags.append(oneflags)

    self.imageflags = imageflags

  # replicate system by Nx and Ny and Nz in each dim
  # rebuild data structs (except imageflags) in this class
  # also alter global boxhi

  def replicate(self,nx,ny,nz):

    natoms = self.natoms
    id = self.id
    label = self.label
    type = self.type
    x = self.x
    y = self.y
    z = self.z
    bonds = self.bonds
    imageflags = self.imageflags

    xbox = boxhi[0]
    ybox = boxhi[1]
    zbox = boxhi[2]

    idnew = []
    labelnew = []
    typenew = []
    xnew = []
    ynew = []
    znew = []
    bondsnew = []

    count = 0
    for k in range(nz):
      for j in range(ny):
        for i in range(nx):
          for m in range(natoms):
            count += 1
            idnew.append(count)
            labelnew.append(label[m])
            typenew.append(type[m])
            xnew.append(str(float(x[m]) + i*xbox))
            ynew.append(str(float(y[m]) + j*ybox))
            znew.append(str(float(z[m]) + k*zbox))

            oneflags = imageflags[m]
            onebonds = []

            for n,mbond in enumerate(bonds[m]):

              # ijk new = which replicated box the mbond atom is in

              if oneflags[n][0] == 0: inew = i
              elif oneflags[n][0] == -1: inew = i - 1
              elif oneflags[n][0] == 1: inew = i + 1
              if inew < 0: inew = nx-1
              if inew == nx: inew = 0

              if oneflags[n][1] == 0: jnew = j
              elif oneflags[n][1] == -1: jnew = j - 1
              elif oneflags[n][1] == 1: jnew = j + 1
              if jnew < 0: jnew = ny-1
              if jnew == ny: jnew = 0

              if oneflags[n][2] == 0: knew = k
              elif oneflags[n][2] == -1: knew = k - 1
              elif oneflags[n][2] == 1: knew = k + 1
              if knew < 0: knew = nz-1
              if knew == nz: knew = 0

              # bondnew = replicated atom ID of bond partner
              # based on replication box (inew,jnew,knew) that owns it

              bondnew = mbond + natoms*(knew*ny*nx + jnew*nx + inew)
              onebonds.append(bondnew)

            bondsnew.append(onebonds)

    self.natoms = natoms * nx*ny*nz
    self.id = idnew
    self.label = labelnew
    self.type = typenew
    self.x = xnew
    self.y = ynew
    self.z = znew
    self.bonds = bondsnew

  # write out new xyzfile for replicated system

  def output(self,outfile):
    fp = open(outfile,'w')
    words = self.header.split()
    print(self.natoms,"replicated",' '.join(words[1:]), file=fp)

    id = self.id
    label = self.label
    type = self.type
    x = self.x
    y = self.y
    z = self.z
    bonds = self.bonds

    # NOTE: worry about formatting of line

    for i in range(self.natoms):
      print(i+1,label[i],x[i],y[i],z[i],type[i], end=' ', file=fp)
      for j in bonds[i]: print(j, end=' ', file=fp)
      print(file=fp)

    fp.close()

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

class PRMfile(object):
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
    self.opbendparams = self.opbend()
    self.ureyparams = self.ureybonds()
    self.pitorsionparams = self.pitorsions()
    self.bitorsionparams = self.bitorsions()
    self.ntypes = len(self.masses)

  # find a section in the PRM file

  def find_section(self,txt):
    txt = "##  %s  ##" % txt
    for iline,line in enumerate(self.lines):
      if txt in line: return iline
    return -1

  # scalar params

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

  # atom classes and masses

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

  # replicate per-atom data: classes and masses

  def replicate_peratom(self,nx,ny,nz):
    classes = self.classes
    masses = self.masses
    natoms = len(self.masses)

    classes_new = []
    masses_new = []

    for k in range(nz):
      for j in range(ny):
        for i in range(nx):
          for m in range(natoms):
            classes_new.append(classes[m])
            masses_new.append(masses[m])

    self.classes = classes_new
    self.masses = masses_new

  # polarization groups

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
    ubflag = 0
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

          option1 = (pflag,ubflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)

          if len(words) >= 7:
            value3 = float(words[6])
            lmp1 = value3
            option2 = (pflag,ubflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)

          if len(words) == 8:
            value4 = float(words[7])
            lmp1 = value4
            option3 = (pflag,ubflag,lmp1,lmp2,lmp3,lmp4,lmp5,lmp6)

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
            # NOTE: just for debugging
            #error("1st bond in BondAngle term not found: %d %d %d" % \
            #      (class1,class2,class3))
            lmp3 = 0.0

          if (class2,class3) in bdict:
            lmp4 = bdict[(class2,class3)]
          elif (class3,class2) in bdict:
            lmp4 = bdict[(class3,class2)]
          else:
            # NOTE: just for debugging
            #error("2nd bond in BondAngle term not found: %d %d %d" % \
            #      (class1,class2,class3))
            lmp4 = 0.0

          params.append((class1,class2,class3,lmp1,lmp2,lmp3,lmp4))
      iline += 1
    return params

  # dihedral interactions

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
            error("torsion has no params: %d %d %d %d" % \
                  (class1,class2,class3,class4))
          if (len(words)-5) % 3:
            error("torsion does not have triplets of params: %d %d %d %d" % \
                  (class1,class2,class3,class4))

          mfourier = int((len(words)-5)/3)
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

  # improper or out-of-plane bend interactions

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

  # convert PRMfile params to LAMMPS angle_style amoeba UB params
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

          params.append((class1,class2,class3,lmp1,lmp2))
      iline += 1
    return params

  # PiTorsion params, will be read from data file by fix pitorsion

  def pitorsions(self):
    params = []
    iline = self.find_section("Pi-Torsion Parameters")
    if iline < 0: return params
    iline += 3

    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "pitors":
          class1 = int(words[1])
          class2 = int(words[2])
          value1 = float(words[3])
          lmp1 = value1

          params.append((class1,class2,lmp1))
      iline += 1
    return params

  # BiTorsion params, will be read from data file by fix bitorsion

  def bitorsions(self):
    params = []
    iline = self.find_section("Torsion-Torsion Parameters")
    if iline < 0: return params
    iline += 3

    while iline < self.nlines:
      words = self.lines[iline].split()
      if len(words):
        if words[0].startswith("###########"): break
        if words[0] == "tortors":
          class1 = int(words[1])
          class2 = int(words[2])
          class3 = int(words[3])
          class4 = int(words[4])
          class5 = int(words[5])
          nx = int(words[6])
          ny = int(words[7])
          iline += 1
          array = []
          for iy in range(ny):
            xrow = []
            for ix in range(nx):
              words = self.lines[iline].split()
              xgrid = float(words[0])
              ygrid = float(words[1])
              value = float(words[2])
              tuple3 = (xgrid,ygrid,value)
              xrow.append(tuple3)
              iline += 1
            array.append(xrow)
          params.append((class1,class2,class3,class4,class5,nx,ny,array))
      iline += 1
    return params

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
bitorsionfile = ""
pbcflag = 0
repflag = 0

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
  elif args[iarg] == "-bitorsion":
    if iarg + 2 > narg: error()
    bitorsionfile = args[iarg+1]
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
  elif args[iarg] == "-rep":
    if iarg + 5 > narg: error()
    repflag = 1
    nxrep = int(args[iarg+1])
    nyrep = int(args[iarg+2])
    nzrep = int(args[iarg+3])
    outxyz = args[iarg+4]
    if outxyz == "NULL": outxyz = ""
    iarg += 5
  else: error()

# error check

if not xyzfile: error("xyzfile not specified")
if not prmfile: error("prmfile not specified")
if not datafile: error("datafile not specified")
if not pbcflag and repflag: error("cannot replicate non-periodic system")
if repflag and (nxrep <= 0 or nyrep <= 0 or nzrep <= 0):
  error("replication factors <= 0 not allowed")

# read Tinker xyz and prm files

xyz = XYZfile(xyzfile)
prm = PRMfile(prmfile)

# replicate system if requested
# both in XYZ and PRM classes

if repflag:
  xyz.bond_images()
  xyz.replicate(nxrep,nyrep,nzrep)
  if outxyz: xyz.output(outxyz)
  prm.replicate_peratom(nxrep,nyrep,nzrep)
  boxhi[0] *= nxrep
  boxhi[1] *= nyrep
  boxhi[2] *= nzrep

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
  for i in range(natoms):
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

# NOTE: could just avoid double count by "if atom1 < atom4" as in bond, angle ?
#       gives different list, but is it still identical ?
#       would have to check 2 data files, write comparison Py script

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

# create pitorsionlist = list of 6-body interactions
# based on central bond, each bond atom is bonded to exactly 2 other atoms
# avoid double counting by requiring atom1 < atom2
# NOTE: need more info on how to order the 6 atoms for Tinker to compute on

type = xyz.type
classes = prm.classes
bonds = xyz.bonds

pitorsionlist = []

for atom1 in id:
  for atom2 in bonds[atom1-1]:
    if atom1 < atom2:
      if len(bonds[atom1-1]) != 3: continue
      if len(bonds[atom2-1]) != 3: continue

      if bonds[atom1-1][0] == atom2:
        atom3 = bonds[atom1-1][1]
        atom4 = bonds[atom1-1][2]
      elif bonds[atom1-1][1] == atom2:
        atom3 = bonds[atom1-1][0]
        atom4 = bonds[atom1-1][2]
      elif bonds[atom1-1][2] == atom2:
        atom3 = bonds[atom1-1][0]
        atom4 = bonds[atom1-1][1]

      if bonds[atom2-1][0] == atom1:
        atom5 = bonds[atom2-1][1]
        atom6 = bonds[atom2-1][2]
      elif bonds[atom2-1][1] == atom1:
        atom5 = bonds[atom2-1][0]
        atom6 = bonds[atom2-1][2]
      elif bonds[atom2-1][2] == atom1:
        atom5 = bonds[atom2-1][0]
        atom6 = bonds[atom2-1][1]

      pitorsionlist.append((atom3,atom4,atom1,atom2,atom5,atom6))

# create bitorsionlist = list of 5-body interactions
# generate topology via double loop over neighbors of central atom3
#   additional double loop over bonds of atom2 and bonds of atom4
# avoid double counting the reverse bitorsion by use of btdict dictionary

type = xyz.type
classes = prm.classes
bonds = xyz.bonds

bitorsionlist = []
btdict = {}

for atom3 in id:
  for atom2 in bonds[atom3-1]:
    for atom4 in bonds[atom3-1]:
      if atom2 == atom4: continue
      for atom1 in bonds[atom2-1]:
        for atom5 in bonds[atom4-1]:
          if atom1 == atom3 or atom5 == atom3 or atom1 == atom5: continue
          if (atom5,atom4,atom3,atom2,atom1) in btdict: continue
          bitorsionlist.append((atom1,atom2,atom3,atom4,atom5))
          btdict[(atom1,atom2,atom3,atom4,atom5)] = 1

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
  else: error("bond not found: %d %d: %d %d" % (atom1,atom2,c1,c2))

  if not flags[m]:
    v1,v2,v3,v4 = params[2:]
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
    #   no need to flip if c1 = c3
    # necessary b/c BondAngle coeffs will be generated with r1,r2 params
    #   from Bond Stretching section of PRM file
    # since in general r1 != r2, the LAMMPS AngleAmoeba class requires
    #   the 3 atoms in the angle be in the order that matches r1 and r2

    if c1 != c3 and (c3,c2,c1) in adict:
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

      #if nbonds != 3:
        #print("Center angle atom has wrong bond count")
        #print("  angle atom IDs:",atom1,atom2,atom3)
        #print("  angle atom classes:",c1,c2,c3)
        #print("  Tinker FF file param options:",len(params[3]))
        #print("  Nbonds and hydrogen count:",nbonds,hcount)
        # NOTE: allow this for now
        #sys.exit()

      if hcount == 0: which = 1
      elif hcount == 1:
        which = 2
        m += 1

      #print("3-bond angle")
      #print("  angle atom IDs:",atom1,atom2,atom3)
      #print("  angle atom classes:",c1,c2,c3)
      #print("  Tinker FF file param options:",len(params[3]))
      #print("  Nbonds and hydrogen count:",nbonds,hcount)
      #print("  which:",which,m)

    elif len(params[3]) == 3:
      nbonds,hcount = xyz.angle_hbond_count(atom1,atom2,atom3,lmptype,lmpmass)

      #if nbonds != 4:
        #print("Center angle atom has wrong bond count")
        #print("  angle atom IDs:",atom1,atom2,atom3)
        #print("  angle atom classes:",c1,c2,c3)
        #print("  Tinker FF file param options:",len(params[3]))
        #print("  Nbonds and hydrogen count:",nbonds,hcount)
        # NOTE: allow this for now
        #sys.exit()

      if hcount == 0: which = 1
      elif hcount == 1:
        which = 2
        m += 1
      elif hcount == 2:
        which = 3
        m += 2

  else:
    error("angle not found: %d %d %d: %d %d %d" % (atom1,atom2,atom3,c1,c2,c3))

  if not flags[m]:
    pflag,ubflag,v1,v2,v3,v4,v5,v6 = params[3][which-1]
    aparams.append((pflag,ubflag,v1,v2,v3,v4,v5,v6))
    flags[m] = len(aparams)
  atype.append(flags[m])

# augment the aparams with bond-angle cross terms from bondangleparams
# generate baparams = LAMMPS bond-angle params for each angle type
# badict = dictionary for angle tuples in bongangleparams

badict = {}
for v1,v2,v3,v4,v5,v6,v7 in prm.bondangleparams:
  if (v1,v2,v3) in badict: continue
  badict[(v1,v2,v3)] = (v4,v5,v6,v7)

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

  if (c1,c2,c3) in badict:
    n1,n2,r1,r2 = badict[(c1,c2,c3)]
  elif (c3,c2,c1) in badict:
    n1,n2,r1,r2 = badict[(c3,c2,c1)]
  else:
    # NOTE: just for debugging
    #print("Bond-stretch angle triplet not found: %d %d %d" % (c1,c2,c3))
    n1,n2,r1,r2 = 4*[0.0]

  baparams.append((n1,n2,r1,r2))

# augment the aparams with Urey_Bradley terms from ureyparams
# generate ubparams = LAMMPS UB params for 1-3 bond in each angle type
# ubdict = dictionary for angle tuples in ureyparams

ubdict = {}
for v1,v2,v3,v4,v5 in prm.ureyparams:
  if (v1,v2,v3) in ubdict: continue
  ubdict[(v1,v2,v3)] = (v4,v5)

ubparams = []

for itype in range(len(aparams)):
  iangle = atype.index(itype+1)
  atom1,atom2,atom3 = alist[iangle]
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]

  # if UB settings exist for this angle type, set ubflag in aparams to 1

  if (c1,c2,c3) in ubdict:
    r1,r2 = ubdict[(c1,c2,c3)]
    pflag,ubflag,v1,v2,v3,v4,v5,v6 = aparams[itype]
    ubflag = 1
    aparams[itype] = (pflag,ubflag,v1,v2,v3,v4,v5,v6)
  elif (c3,c2,c1) in ubdict:
    r1,r2 = ubdict[(c3,c2,c1)]
    pflag,ubflag,v1,v2,v3,v4,v5,v6 = aparams[itype]
    ubflag = 1
    aparams[itype] = (pflag,ubflag,v1,v2,v3,v4,v5,v6)
  else:
    r1,r2 = 2*[0.0]

  ubparams.append((r1,r2))

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
    error("dihedral not found: %d %d %d %d: %d %d %d %d" % \
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

# generate pitorsiontype = LAMMPS type of each pitorsion
# generate pitorsionparams = LAMMPS params for each pitorsion type
# flags[i] = which LAMMPS pitorsion type (1-N) the Ith Tinker PRM file ptors is
#        0 = none
# convert prm.pitorsionparams to a dictionary for efficient searching
# key = (class1,class2)
# value = (M,params) where M is index into prm.pitorsionparams

id = xyz.id
type = xyz.type
classes = prm.classes

pitdict = {}
for m,params in enumerate(prm.pitorsionparams):
  pitdict[(params[0],params[1])] = (m,params)

flags = len(prm.pitorsionparams)*[0]
pitorsiontype = []
pitorsionparams = []
pitorsionlist_reduced = []

for tmp1,tmp2,atom1,atom2,tmp3,tmp4 in pitorsionlist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]

  # 6-tuple is only a PiTorsion if central 2 atoms match an entry in PRM file
  # pitorsionlist_reduced = list of just these 6-tuples

  if (c1,c2) in pitdict or (c2,c1) in pitdict:
    if (c1,c2) in pitdict: m,params = pitdict[(c1,c2)]
    else: m,params = pitdict[(c2,c1)]
    pitorsionlist_reduced.append((tmp1,tmp2,atom1,atom2,tmp3,tmp4))

    if not flags[m]:
      v1 = params[2:]
      pitorsionparams.append(v1)
      flags[m] = len(pitorsionparams)
    pitorsiontype.append(flags[m])

# replace original pitorsionlist with reduced version

pitorsionlist = pitorsionlist_reduced

# generate bitorsiontype = LAMMPS type of each bitorsion
# generate bitorsionparams = LAMMPS params for each bitorsion type
# flags[i] = which LAMMPS bitorsion type (1-N) the Ith Tinker PRM file btors is
#        0 = none
# convert prm.bitorsionparams to a dictionary for efficient searching
# key = (class1,class2,class3,class4,class5)
# value = (M,params) where M is index into prm.bitorsionparams

id = xyz.id
type = xyz.type
classes = prm.classes

bitdict = {}
for m,params in enumerate(prm.bitorsionparams):
  bitdict[(params[0],params[1],params[2],params[3],params[4])] = (m,params)

flags = len(prm.bitorsionparams)*[0]
bitorsiontype = []
bitorsionparams = []
bitorsionlist_reduced = []

for atom1,atom2,atom3,atom4,atom5 in bitorsionlist:
  type1 = type[atom1-1]
  type2 = type[atom2-1]
  type3 = type[atom3-1]
  type4 = type[atom4-1]
  type5 = type[atom5-1]
  c1 = classes[type1-1]
  c2 = classes[type2-1]
  c3 = classes[type3-1]
  c4 = classes[type4-1]
  c5 = classes[type5-1]

  # 5-tuple is only a BiTorsion if 5 atoms match an entry in PRM file
  # bitorsionlist_reduced = list of just these 5-tuples

  if (c1,c2,c3,c4,c5) in bitdict or (c5,c4,c3,c2,c1) in bitdict:
    if (c1,c2,c3,c4,c5) in bitdict: m,params = bitdict[(c1,c2,c3,c4,c5)]
    else: m,params = bitdict[(c5,c4,c3,c2,c1)]
    bitorsionlist_reduced.append((atom1,atom2,atom3,atom4,atom5))

    if not flags[m]:
      v1 = params[5:]
      bitorsionparams.append(v1)
      flags[m] = len(bitorsionparams)
    bitorsiontype.append(flags[m])

# replace original bitorsionlist with reduced version

bitorsionlist = bitorsionlist_reduced

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
npitorsions = len(pitorsionlist)
nbitorsions = len(bitorsionlist)

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
  for i,one in enumerate(ubparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["UreyBradley Coeffs"] = lines

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

if npitorsions:
  d.headers["pitorsions"] = len(pitorsionlist)
  d.headers["pitorsion types"] = len(pitorsionparams)

  lines = []
  for i,one in enumerate(pitorsionparams):
    strone = [str(single) for single in one]
    line = "%d %s" % (i+1,' '.join(strone))
    lines.append(line+'\n')
  d.sections["PiTorsion Coeffs"] = lines

  lines = []
  for i,one in enumerate(pitorsionlist):
    line = "%d %d %d %d %d %d %d %d" % \
           (i+1,pitorsiontype[i],one[0],one[1],one[2],one[3],one[4],one[5])
    lines.append(line+'\n')
  d.sections["PiTorsions"] = lines

if nbitorsions:
  d.headers["bitorsions"] = len(bitorsionlist)

  # if there are bitorsions, then -bitorsion file must have been specified

  if not bitorsionfile:
    error("no -bitorsion file was specified, but %d bitorsions exist" % \
          nbitorsions)

  fp = open(bitorsionfile,'w')
  print("Tinker BiTorsion parameter file for fix bitorsion\n", file=fp)
  print("%d bitorsion types" % len(bitorsionparams), file=fp)
  itype = 0
  for nx,ny,array in bitorsionparams:
    itype += 1
    print(file=fp)
    print(itype,nx,ny, file=fp)
    for ix in range(nx):
      for iy in range(ny):
        xgrid,ygrid,value = array[ix][iy]
        print(" ",xgrid,ygrid,value, file=fp)
  fp.close()

  lines = []
  for i,one in enumerate(bitorsionlist):
    line = "%d %d %d %d %d %d %d" % \
           (i+1,bitorsiontype[i],one[0],one[1],one[2],one[3],one[4])
    lines.append(line+'\n')
  d.sections["BiTorsions"] = lines

d.write(datafile)

# print stats to screen

print("Natoms =",natoms)
print("Ntypes =",ntypes)
print("Tinker XYZ types =",len(tink2lmp))
print("Tinker PRM types =",prm.ntypes)
#print("Tinker groups =",ngroups)
print("Nmol =",nmol)
print("Nbonds =",nbonds)
print("Nangles =",nangles)
print("Ndihedrals =",ndihedrals)
print("Nimpropers =",nimpropers)
print("Npitorsions =",npitorsions)
print("Nbitorsions =",nbitorsions)
print("Nbondtypes =",len(bparams))
print("Nangletypes =",len(aparams))
print("Ndihedraltypes =",len(dparams))
print("Nimpropertypes =",len(oparams))
print("Npitorsiontypes =",len(pitorsionparams))
print("Nbitorsiontypes =",len(bitorsionparams))
