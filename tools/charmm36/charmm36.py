
# python3 charmm36.py | gzip -9 > charmm_c36_jul24.gz

################################################################
# FIXME: dihedral weighting factor
################################################################

import re

mass = dict()

# BONDS
# V(bond) = Kb(b - b0)**2
# Kb: kcal/mole/A**2
# b0: A
bond = dict()

# ANGLES
# V(angle) = Ktheta(Theta - Theta0)**2
# V(Urey-Bradley) = Kub(S - S0)**2
# Ktheta: kcal/mole/rad**2
# Theta0: degrees
# Kub: kcal/mole/A**2 (Urey-Bradley)
# S0: A
angle = dict()

# DIHEDRALS
# V(dihedral) = Kchi(1 + cos(n(chi) - delta))
# Kchi: kcal/mole
# n: multiplicity
# delta: degrees
dihedral = dict()

# IMPROPER
# V(improper) = Kpsi(psi - psi0)**2
# Kpsi: kcal/mole/rad**2
# psi0: degrees
# note that the second column of numbers (0) is ignored
improper = dict()

# NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
# cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
# V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
# epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
# Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
# atom  ignored    epsilon      Rmin/2   ignored   eps,1-4       Rmin/2,1-4
pair = dict()

prms = [
  "par_all36m_prot.prm",
  "par_all36_na.prm",
  #"par_all36_carb.prm",
  "par_all36_lipid.prm",
  "par_all36_cgenff.prm",
  #"toppar_all36_moreions.str",
  #"toppar/par_interface.prm",
  "toppar_water_ions.str"]

#prms = ["par_all36_lipid.prm"]

for prm in prms:

  file = open(prm, "r")

  for line in file:

    match = re.search(r"^MASS\s+-1\s+(\w+)\s+(-?\d+.\d+).*", line)
    if( match != None ):
      mass.update( {match.group(1): match.group(2)} )

    match = re.search(r"^(\w+)\s+(\w+)\s+(\d+.\d+)\s+(\d+.\d+)\s+.*", line)
    if( match != None ):
      bond.update( {"{}-{}".format(match.group(1),match.group(2)) :
      "{} {}".format(match.group(3),match.group(4))} )

    match = re.search(r"^(\w+)\s+(\w+)\s+(\w+)\s+(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+)\s+(\d+.\d+).*", line)
    if( match != None ):
      angle.update( {"{}-{}-{}".format(match.group(1),match.group(2),match.group(3)) :
      "{} {} {} {}".format(match.group(4),match.group(5),match.group(6),match.group(7))} )

    match = re.search(r"^(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(-?\d+.\d+)\s+(\d+)\s+(\d+).*", line)
    if( match != None ):
      dihedral.update( {"{}-{}-{}-{}".format(match.group(1),match.group(2),match.group(3),match.group(4)) :
      "{} {} {} 1.00".format(match.group(5),match.group(6),match.group(7))} )

    match = re.search(r"^(\w+)\s+(\w+)\s+(\w+)\s+(\w+)\s+(\d+.\d+)\s+0\s+(\d+.\d+).*", line)
    if( match != None ):
      improper.update( {"{}-{}-{}-{}".format(match.group(1),match.group(2),match.group(3),match.group(4)) :
      "{} {}".format(match.group(5),match.group(6))} )

  file.close()


#44 atoms
#11 atom types
#42 bonds
#15 bond types
#74 angles
#29 angle types
#100 dihedrals
#36 dihedral types
#44 impropers
#13 improper types

# Header

print( "LAMMPS CHARMM36 force field (toppar_c36_jul24.tgz) [https://mackerell.umaryland.edu/charmm_ff.shtml]\n" )

print( "  ", len(mass), " atom types" )
print( "  ", len(bond), " bond types" )
print( "  ", len(angle), " angle types" )
print( "  ", len(dihedral), " dihedral types" )
print( "  ", len(improper), " improper types" )

#  -------- Atom Type Labels  --------
print( "\nAtom Type Labels\n" )
i=1
for k in mass.keys():
  print("  ", i, k)
  i+=1

# -------- Masses --------
print( "\nMasses\n" )
i=1
for v  in mass.values():
  print("  ", i, v)
  i+=1

# -------- Bond Type Labels --------
print( "\nBond Type Labels\n" )
i=1
for k in bond.keys():
  print("  ", i, k)
  i+=1

# -------- Bond Coeffs --------
print( "\nBond Coeffs # harmonic\n" )
i=1
for v in bond.values():
  print("  ", i, v)
  i+=1

# -------- Angle Type Labels --------
print( "\nAngle Type Labels\n" )
i=1
for k in angle.keys():
  print("  ", i, k)
  i+=1

# -------- Angle Coeffs --------
print( "\nAngle Coeffs # charmm\n" )
i=1
for v in angle.values():
  print("  ", i, v)
  i+=1

# -------- Dihedral Type Labels --------
print( "\nDihedral Type Labels\n" )
i=1
for k in dihedral.keys():
  print("  ", i, k)
  i+=1

# -------- Dihedral Coeffs --------
print( "\nDihedral Coeffs # charmmfsw\n" )
i=1
for v in dihedral.values():
  print("  ", i, v)
  i+=1

# -------- Improper Type Labels --------
print( "\nImproper Type Labels\n" )
i=1
for k in improper.keys():
  print("  ", i, k)
  i+=1

# -------- Improper Coeffs --------
print( "\nImproper Coeffs # harmonic\n" )
i=1
for v in improper.values():
  print("  ", i, v)
  i+=1




# -------- Pair Coeffs --------
