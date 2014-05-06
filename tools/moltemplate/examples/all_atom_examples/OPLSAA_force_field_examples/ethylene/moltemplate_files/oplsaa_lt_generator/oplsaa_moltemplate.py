#! /usr/bin/env python
#
# The purpose of this script is to create a moltemplate lt file for the oplsaa.
# forcefield.  This will assist researchers in building complex simulations using
# this OPLS-UA and the OPLS-AA forcefields.

__author__="Jason Lambert"
# (some additional corrections by Miguel Gonzalez and Andrew Jewett)
__version__="0.16"

import sys
import os
from operator import itemgetter

print("""

Warning:
  Run this program on a SUBSET of the OPLS atoms which are relevant
  to your problem.  If not, this program (and moltemplate) may crash 
  your computer and/or generate enormous files that you do not need.
""")
#  To do that, first make a copy of the \"oplsaa.prm\" file
#  (which can be downloaded from the TINKER web site).
#  The lines in this file beginning with the word \"atoms\" should
#  define the atoms which you plan to put in your simulation. All other 
#  lines beginning with the word \"atoms\" should be deleted.
#  (Leave the other sections of this file alone.)
#""")



#input data from file containing opls aa force field parameters.
try:
   f=open(sys.argv[1],"r")
except:
   print("need to specify file name as an input argument:")
   print("python oplsaa_moltemplate.py <forcefield file name>")
   print("or file name is specified incorrectly")
   sys.exit()
#output lt file
g=open("oplsaa.lt","w")



lines = f.readlines()



# Ignore/Comment out lines before the "##  Atom Type Definitions  ##" section.

for i in range(0, len(lines)):
   if (lines[i].find("##  Atom Type Definitions  ##") != -1):
      break
   else:
      lines[i] = '# ' + lines[i]


# As of 2014-4-19, there appear to be 906 atom types, but we don't assume this.
# First try to infer out how many atom types there were in the original 
# oplsaa.prm file, or at least find an upper bound on the atom-type numbers.
# (Keep track of the maximum value of the first column in the "atom" section.)
max_atomType = 0
for line in lines:
   # skip over text after a # comment character
   ic = line.find('#') 
   if ic != -1:
      line = (line[:ic]).strip()
   else:
      line = line.strip()
   # now look for lines beginning with the word "atom"
   tokens = line.split()
   if ((len(tokens)>2) and 
       (tokens[0] == "atom") and 
       (int(tokens[1]) > max_atomType)):
      max_atomType = int(tokens[1])


#temporary storage file
h=open("temp.txt","w+")
atom_lookup={} #this dictionary contains all the atom ffid's as a key and the number of atoms with that key
#atom=[[10000,10000] for i in range(906)]  <- don't assume there are 906 atoms
atom=[[-10000,-10000] for i in range(0,max_atomType+1)]
#charge_by_type={} # lookup charge by atom type
#vdw_by_type={}    # lookup epsilon & sigma paramters by atom type
charge_by_type=[0.0 for i in range(0,max_atomType+1)]    # lookup charge by atom
vdw_by_type=[(0.0,0.0) for i in range(0,max_atomType+1)] # lookup epsilon & sigma



#atom is declared this way so for sorting purposes.  
#atom contains the following data upon allocation
#atom[][0]=atom_id( Important for partial charges and non_bonded interactions)
#atom[][1]=atom_ffid( Important for stretches, bending, torsions and impropers)
#atom[][2]=atom_mass
#atom[][3]=partial charge
#atom[][4]=non_bonding sigma
#atom[][5]=non_bonding epsilon
#atom[][6]=atom comment
bond=[]
#bond contains the following data
#bond[0]=atom 1 ffid
#bond[1]=atom 2 ffid
#bond[2]=bond spring constant(OPLS-aa compatible)
#bond[3]=equilibrium bond distance(Angstrom)
angle=[]
#angle contains the following data
#angle[0]=atom 1 ffid
#angle[1]=atom 2 ffid
#angle[2]=atom 3 ffid
#angle[3]=spring constant
#angle[4]=equilibrium angle (degrees)
dihedral=[]
#dihedral contains the following data
#dihedral[0]=atom 1 ffid
#dihedral[1]=atom 2 ffid
#dihedral[2]=atom 3 ffid
#dihedral[3]=atom 4 ffid
#dihedral[4]=v1
#dihedral[5]=v2
#dihedral[6]=v3
#dihedral[7]=v4
improper=[]
#improper[0]=atom 1 ffid
#improper[1]=atom 2 ffid(central atom)
#improper[2]=atom 3 ffid
#improper[3]=atom 4 ffid
#improper[4]=spring coefficient
#improper[5]=equilibrium angle


#This section gets all the parameters from the force field file
for line in lines:

   # skip over text after a # comment character
   ic = line.find('#') 
   if ic != -1:
      line = (line[:ic]).strip()
   else:
      line = line.strip()

   if line.find("atom") == 0:
      line=line.split()
      atom[int(line[1])-1]=[int(line[1]),int(line[2]),float(line[-2]),
                            0.0,0.0,0.0," ".join(line[3:-2])]    
   elif line.find("vdw") == 0:
      line=line.split()
      #vdw_temp.append([float(line[1]),float(line[2]),float(line[3])])
      if (int(line[1]) <= max_atomType):
         vdw_by_type[int(line[1])] = (float(line[2]),float(line[3]))
   elif line.find("bond") == 0:
      line=line.split()
      bond.append([int(line[1]),int(line[2]),float(line[3]),float(line[4])])
   elif line.find("angle") == 0:
      line=line.split()
      angle.append([int(line[1]),int(line[2]),int(line[3]),
      float(line[4]),float(line[5])])
   elif line.find("torsion") == 0:
      line=line.split()
      dihedral.append([int(line[1]),int(line[2]),int(line[3]),int(line[4]),
      float(line[5]),float(line[8]), float(line[11]), 0.0])   	
   elif line.find("charge") == 0:
      line=line.split()
      #charge_temp.append([int(line[1]),float(line[2])])
      if (int(line[1]) <= max_atomType):
         charge_by_type[int(line[1])] = float(line[2])
   elif line.find("imptors") == 0:
      line=line.split()
      improper.append([int(line[1]), int(line[2]), 
      int(line[3]), int(line[4]), float(line[5]), float(line[6])])

if len(atom) > 600:
   sys.stderr.write("WARNING: The number of atom types in your file exceeds 600\n"
                    "         (You were supposed to edit out the atoms you don't need.\n"
                    "          Not doing this may crash your computer.)\n"
                    "\n"
                    "         Proceed? (Y/N): ")
   reply = sys.stdin.readline()
   if find(reply.strip().lower(), 'y') != 0:
      exit(0)

#adding the charge and Lennard Jones parameters to
#to each atom type.
#----------------------------------------------#
for i in range(0,len(atom)):
   atom_type_num = atom[i][0]
   #q = charge_by_type.get(atomTypeNum)
   #if q:
   #   atom[i][3] = q
   if atom_type_num != -10000:
      q = charge_by_type[atom_type_num]
      atom[i][3] = q

for i in range(0,len(atom)):
   atom_type_num = atom[i][0]
   #vdw_params = vdw_by_type.get(atomTypeNum)
   #if vdw_params:
   #   atom[i][4] = vdw_params[0]
   #   atom[i][5] = vdw_params[1]
   if atom_type_num != -10000:
      vdw_params = vdw_by_type[atom_type_num]
      atom[i][4] = vdw_params[0]
      atom[i][5] = vdw_params[1]
   
del(charge_by_type)
del(vdw_by_type)

#----------------------------------------------------------#
#begin writing content to lt file
g.write("OPLSAA {\n\n" )

#write out the atom masses
#----------------------------------------------------------#
g.write("  write_once(\"Data Masses\"){\n")#checked with gaff
for i,x in enumerate(atom):
   if x[0] != -10000:
     g.write("    @atom:{} {} #{} partial charge={}\n".format(
     x[0],x[2],x[6],x[3]))
g.write("  } #(end of atom masses)\n\n")
#----------------------------------------------------------#


#write out the pair coefficients
#----------------------------------------------------------#
g.write("  write_once(\"In Settings\"){\n")#checked with gaff
for i,x in enumerate(atom):
  if x[0] != -10000:
   g.write("    pair_coeff @atom:{0} @atom:{0} lj/cut/coul/long {1} {2}\n".format(x[0],x[5],x[4]))
g.write("  } #(end of pair coeffs)\n\n")

g.write("  write_once(\"In Charges\"){\n")#checked with gaff
for i,x in enumerate(atom):
  if x[0] != -10000:
   g.write("    set type @atom:{0} charge {1}\n".format(x[0],x[3]))
g.write("  } #(end of atom charges)\n\n")

#-----------------------------------------------------------#

# This part of the code creates a lookup dictionary
# that allows you to find every type of atom by its 
# force field id. force field id is the id number
# relevant to bonds, angles, dihedrals, and impropers.
# This greatly increases the speed of angle, bond, dihedral
# and improper assignment.
#------------------------------------------------------------#
atom=sorted(atom,key=itemgetter(1))
atom_ffid=0
for x in atom:
	if x[1]==atom_ffid:
		atom_lookup[x[1]].append(x[0])
	elif x[1]>atom_ffid:
	   atom_lookup[x[1]]=[x[0]]
	   atom_ffid=x[1]
atom_lookup[0]=["*"]
#-------------------------------------------------------------#


#writing out the bond coefficients and bond parameters#
#-------------------------------------------------------------#	   
g.write("  write_once(\"In Settings\") {\n")
index1=0
for x in bond:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      #g.write("    bond_coeff @bond:{}-{} harmonic {} {}\n".format(y,z,x[2]/2,x[3]))
      # Miguel Gonzales corrected this line to:
      g.write("    bond_coeff @bond:{}-{} harmonic {} {}\n".format(y,z,x[2],x[3]))
      h.write("    @bond:{0}-{1} @atom:{0} @atom:{1}\n".format(y,z))
g.write("  } #(end of bond_coeffs)\n\n")
h.seek(0,0)      
g.write("  write_once(\"Data Bonds By Type\") {\n")
for line in h.readlines():
  g.write(line)
g.write("  } #(end of bonds by type)\n\n")
del(bond)
h.close()
#-----------------------------------------------------------#
h=open("temp.txt","w+")

#writing out angle coefficients and angles by type.---------#
#-----------------------------------------------------------#  
g.write("  write_once(\"Data Angles By Type\"){\n") 
for x in angle:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]):
         #print(y,z,u,x)
         #h.write("    angle_coeff @angle:{}-{}-{} harmonic {} {}\n".format(y,z,u,x[3]/2.0,x[4]))
         # Miguel Gonzales corrected this line:
         h.write("    angle_coeff @angle:{}-{}-{} harmonic {} {}\n".format(y,z,u,x[3],x[4]))
         g.write("    @angle:{0}-{1}-{2} @atom:{0} @atom:{1} @atom:{2}\n".format(y,z,u))


g.write("  } #(end of angles by type)\n\n") 
h.seek(0,0)
g.write("  write_once(\"In Settings\" ){\n")
for line in h.readlines():
	g.write(line)
g.write("  } #(end of angle_coeffs)\n\n")                    
del(angle)           
h.close()
#----------------------------------------------------------#

#writing dihedrals by type and dihedral coefficients-------#
h=h=open("temp.txt","w+")
g.write("  write_once(\"Data Dihedrals By Type\") {\n")
#print(atom_lookup)
for x in dihedral:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]): 
        for v in atom_lookup.get(x[3],[]):
          if x[0]!=0 and x[3]!=0:
            g.write("    @dihedral:{0}-{1}-{2}-{3} @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            h.write("    dihedral_coeff @dihedral:{}-{}-{}-{} opls {} {} {} {}\n".format(
            y,z,u,v,x[4],x[5],x[6],x[7]))
          elif x[0]==0 and x[3]!=0:
            g.write("    @dihedral:0-{1}-{2}-{3} @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            h.write("    dihedral_coeff @dihedral:0-{}-{}-{} opls {} {} {} {}\n".format(
            z,u,v,x[4],x[5],x[6],x[7]))
          elif x[0]==0 and x[3]==0:
            g.write("    @dihedral:0-{1}-{2}-0 @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            #h.write("    dihedral_coeff @dihedral:0-{}-{}-0 harmonic {} {} {} {}\n".format(
            h.write("    dihedral_coeff @dihedral:0-{}-{}-0 opls {} {} {} {}\n".format(
            z,u,x[4],x[5],x[6],x[7]))       
 
del(dihedral)
g.write("  } #(end of Dihedrals by type)\n\n")               
h.seek(0,0)   
g.write("  write_once(\"In Settings\") {\n") 
for line in h.readlines():
   g.write(line)
g.write("  } #(end of dihedral_coeffs)\n\n")
h.close()
#-----------------------------------------------------------------------#

#----writing out improper coefficients and impropers by type------------# 
h=open("temp.txt","w+")
g.write("  write_once(\"Data Impropers By Type\") {\n")
for x in improper:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]): 
        for v in atom_lookup.get(x[3],[]):
         if x[0]==0 and x[1]==0 and x[3]==0:
            g.write("    @improper:{2}-0-0-0 @atom:{2} @atom:{0} @atom:{1} @atom:{3}\n".format(
            y,z,u,v))
            h.write("    improper_coeff @improper:{2}-0-0-0 harmonic {4} {5} \n".format(
            y,z,u,v,x[4]/2,0))
         else:
            g.write("    @improper:{2}-0-0-{3} @atom:{2} @atom:{0} @atom:{1} @atom:{3}\n".format(
            y,z,u,v))
            h.write("    improper_coeff @improper:{2}-0-0-{3} harmonic {4} {5} \n".format(
            y,z,u,v,x[4]/2,0))
            
                
g.write("  } #(end of impropers by type)\n\n") 
h.seek(0,0)   
g.write("  write_once(\"In Settings\") {\n") 
for line in h.readlines():
   g.write(line)
g.write("  } #(end of improp_coeffs)\n\n")
#-----------------------------------------------------------------------#

#This section writes out the input parameters required for an opls-aa simulation
# lammps.
g.write("  write_once(\"In Init\") {\n")
g.write("    units real\n")
g.write("    atom_style full\n")
g.write("    bond_style hybrid harmonic\n")
g.write("    angle_style hybrid harmonic\n")
g.write("    dihedral_style hybrid opls\n")
g.write("    improper_style hybrid harmonic\n")
#g.write("    pair_style hybrid lj/cut/coul/cut 10.0 10.0\n")
g.write("    pair_style hybrid lj/cut/coul/long 10.0 10.0\n")
g.write("    pair_modify mix arithmetic\n")
g.write("    special_bonds lj/coul 0.0 0.0 0.5\n")    
g.write("    kspace_style pppm 0.0001\n")
g.write("  } #end of init parameters\n\n")
g.write("} # OPLSAA\n")
f.close()
g.close()
h.close()
os.remove("temp.txt")





   
