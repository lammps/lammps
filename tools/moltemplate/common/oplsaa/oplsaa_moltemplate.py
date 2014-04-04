#! /usr/bin/env python
#The purpose of this script is to create a moltemplate lt file for the opls-aa forcefield.
#This will assist researchers in building complex simulations using this OPLS-UA and the OPLS-AA forcefields.
__author__="Jason Lambert"
__version__="0.13"

import sys
import os
from operator import itemgetter

print("""

Warning:
  Run this program on a SUBSET of the OPLS forcefield relevant
  to your problem.  It is possible for you to generate a full OPLS force
  field moltemplate file, but that demands a lot of time and generates
  a file that is on the order 100 gigabytes (which is too large for
  moltemplate to read).
  Save yourself time and energy, make a copy of the oplsaa.txt that 
  only contains atoms (in the \"atoms\" section) relevant to your problem.    
""")



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

#temporary storage file
h=open("temp.txt","w+")
atom_lookup={} #this dictionary contains all the atom ffid's as a key and the number of atoms with that key
atom=[[10000,10000] for i in range(906)]
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

charge_temp=[] #temporarily store the charges
vdw_temp=[]


#This section gets all the parameters from the force field file
for line in f.readlines():
   if "atom" in line and "#" not in line:
	   line=line.split()
	   atom[int(line[1])-1]=[int(line[1]),int(line[2]),float(line[-2]),
	   0.0,0.0,0.0," ".join(line[3:-2])]      
   elif "vdw" in line and "#" not in line:
      line=line.split()
      vdw_temp.append([float(line[1]),float(line[2]),float(line[3])])
   elif "bond" in line and "#" not in line:
      line=line.split()
      bond.append([int(line[1]),int(line[2]),float(line[3]),float(line[4])])
   elif "angle" in line and "#" not in line:
      line=line.split()
      angle.append([int(line[1]),int(line[2]),int(line[3]),
      float(line[4]),float(line[5])])
   elif "torsion" in line and "#" not in line:
      line=line.split()
      dihedral.append([int(line[1]),int(line[2]),int(line[3]),int(line[4]),
      float(line[5]),float(line[8]), float(line[11]), 0.0])   	
   elif "charge" in line and "#" not in line:
      line=line.split()
      charge_temp.append([int(line[1]),float(line[2])])
   elif "imptors" in line and "#" not in line:
      line=line.split()
      improper.append([int(line[1]), int(line[2]), 
      int(line[3]), int(line[4]), float(line[5]), float(line[6])])

#adding the charge and Lennard Jones parameters to
#to each atom type.
#----------------------------------------------#
i=0
for j,x in enumerate(charge_temp):
  if x[0]==atom[i][0]:
    atom[i][3]=x[1]
  i=i+1
    #print(x[1])
i=0
for j,x in enumerate(vdw_temp):
  #print x
  if x[0]==atom[i][0]:
    atom[i][4]=x[1]
    atom[i][5]=x[2]
  i=i+1
del(charge_temp)
del(vdw_temp)
#----------------------------------------------------------#
#begin writing content to lt file
g.write("OPLSAA {\n\n\n" )

#write out the atom masses
#----------------------------------------------------------#
g.write("write_once(\"Data Masses\"){\n")#checked with gaff
for i,x in enumerate(atom):
   if x[0]<10000:              
     g.write("@atom:{} {} #{} partial charge={}\n".format(
     x[0],x[2],x[6],x[3]))
g.write("} #(end of atom masses)\n\n")
#----------------------------------------------------------#


#write out the pair coefficients
#----------------------------------------------------------#
g.write("write_once(\"In Settings\"){\n")#checked with gaff
for i,x in enumerate(atom):
  if x[0]<10000:
   g.write("pair_coeff @atom:{0} @atom:{0} lj/cut/coul/long {1} {2}\n".format(x[0],x[5],x[4]))
g.write("} #(end of pair coeffs)\n\n")

g.write("write_once(\"In Charges\"){\n")#checked with gaff
for i,x in enumerate(atom):
  if x[0]<10000:
   g.write("set type @atom:{0} charge {1}\n".format(x[0],x[3]))
g.write("} #(end of atom charges)\n\n")

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
g.write("write_once(\"In Settings\") {\n")
index1=0
for x in bond:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      g.write("bond_coeff @bond:{}-{} harmonic {} {}\n".format(y,z,x[2]/2,x[3]))
      h.write("@bond:{0}-{1} @atom:{0} @atom:{1}\n".format(y,z))
g.write("} #(end of bond_coeffs)\n\n")
h.seek(0,0)      
g.write("write_once(\"Data Bonds By Type\") {\n")
for line in h.readlines():
  g.write(line)
g.write("} #(end of bonds by type)\n\n")
del(bond)
h.close()
#-----------------------------------------------------------#
h=open("temp.txt","w+")

#writing out angle coefficients and angles by type.---------#
#-----------------------------------------------------------#  
g.write("write_once(\"Data Angles By Type\"){\n") 
for x in angle:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]):
         #print(y,z,u,x)
            h.write("angle_coeff @angle:{}-{}-{} harmonic {} {}\n".format(y,z,u,
            x[3]/2.0,x[4]))
            g.write("@angle:{0}-{1}-{2} @atom:{0} @atom:{1} @atom:{2}\n".format(
            y,z,u))


g.write("} #(end of angles by type)\n\n") 
h.seek(0,0)
g.write("write_once(\"In Settings\" ){\n")
for line in h.readlines():
	g.write(line)
g.write("} #(end of angle_coeffs)\n\n")                    
del(angle)           
h.close()
#----------------------------------------------------------#

#writing dihedrals by type and dihedral coefficients-------#
h=h=open("temp.txt","w+")
g.write("write_once(\"Data Dihedrals By Type\") {\n")
#print(atom_lookup)
for x in dihedral:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]): 
        for v in atom_lookup.get(x[3],[]):
          if x[0]!=0 and x[3]!=0:
            g.write("@dihedral:{0}-{1}-{2}-{3} @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            h.write("dihedral_coeff @dihedral:{}-{}-{}-{} opls {} {} {} {}\n".format(
            y,z,u,v,x[4],x[5],x[6],x[7]))
          elif x[0]==0 and x[3]!=0:
            g.write("@dihedral:0-{1}-{2}-{3} @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            h.write("dihedral_coeff @dihedral:0-{}-{}-{} opls {} {} {} {}\n".format(
            z,u,v,x[4],x[5],x[6],x[7]))
          elif x[0]==0 and x[3]==0:
            g.write("@dihedral:0-{1}-{2}-0 @atom:{0} @atom:{1} @atom:{2} @atom:{3}\n".format(
            y,z,u,v))
            #h.write("dihedral_coeff @dihedral:0-{}-{}-0 harmonic {} {} {} {}\n".format(
            h.write("dihedral_coeff @dihedral:0-{}-{}-0 opls {} {} {} {}\n".format(
            z,u,x[4],x[5],x[6],x[7]))       
 
del(dihedral)
g.write("} #(end of Dihedrals by type)\n\n")               
h.seek(0,0)   
g.write("write_once(\"In Settings\") {\n") 
for line in h.readlines():
   g.write(line)
g.write("} #(end of dihedral_coeffs)\n\n")
h.close()
#-----------------------------------------------------------------------#

#----writing out improper coefficients and impropers by type------------# 
h=open("temp.txt","w+")
g.write("write_once(\"Data Impropers By Type\") {\n")
for x in improper:
  for y in atom_lookup.get(x[0],[]):
    for z in atom_lookup.get(x[1],[]):
      for u in atom_lookup.get(x[2],[]): 
        for v in atom_lookup.get(x[3],[]):
         if x[0]==0 and x[1]==0 and x[3]==0:
            g.write("@improper:{2}-0-0-0 @atom:{2} @atom:{0} @atom:{1} @atom:{3}\n".format(
            y,z,u,v))
            h.write("improper_coeff @improper:{2}-0-0-0 harmonic {4} {5} \n".format(
            y,z,u,v,x[4]/2,0))
         else:
            g.write("@improper:{2}-0-0-{3} @atom:{2} @atom:{0} @atom:{1} @atom:{3}\n".format(
            y,z,u,v))
            h.write("improper_coeff @improper:{2}-0-0-{3} harmonic {4} {5} \n".format(
            y,z,u,v,x[4]/2,0))
            
                
g.write("} #(end of impropers by type)\n\n") 
h.seek(0,0)   
g.write("write_once(\"In Settings\") {\n") 
for line in h.readlines():
   g.write(line)
g.write("} #(end of improp_coeffs)\n\n")
#-----------------------------------------------------------------------#

#This section writes out the input parameters required for an opls-aa simulation
# lammps.
g.write("write_once(\"In Init\") {\n")
g.write("units real\n")
g.write("atom_style full\n")
g.write("bond_style hybrid harmonic\n")
g.write("angle_style hybrid harmonic\n")
g.write("dihedral_style hybrid opls\n")
g.write("improper_style hybrid harmonic\n")
#g.write("pair_style hybrid lj/cut/coul/cut 10.0 10.0\n")
g.write("pair_style hybrid lj/cut/coul/long 10.0 10.0\n")
g.write("pair_modify mix arithmetic\n")
g.write("special_bonds lj 0.0 0.0 0.5\n")    
g.write("kspace_style pppm 0.0001\n")
g.write("} #end of init parameters\n")
g.write("} # OPLSAA\n")
f.close()
g.close()
h.close()
os.remove("temp.txt")



      

   
   




   
