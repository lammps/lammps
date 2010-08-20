Info="""
Module name: lmp2xyz.py 

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@caltech.edu
Project: pEFF
Version: August 2009

Extracts the xyz from a lammps trajectory dump of style custom:
dump    1 all custom period dump_file id type x y z spin radius ...

Usage: python lmp2xyz.py lammps_dump_filename xyz_filename

"""

import os, sys
from math import log10
from numpy import zeros

mass={"1.00794":"H","4.002602":"He","6.941":"Li","9.012182":"Be","10.811":"B","12.0107":"C","1.00":"Au","0.0005486":"Au"}

def lmp2xyz(lammps,xyz):
  print "\nGenerating %s file"%(xyz)
  fin=open(lammps,'r')
  fout=open(xyz,'w')
  header=9
  lines=fin.readlines()
  numatoms=lines[3].split()[0]
  fsize=os.system("wc -l %s> lines"%(lammps))
  tmp=open('lines','r')
  tlines=tmp.readline()
  tmp.close()
  flines=int(tlines.split()[0])
  snaps=flines/(int(numatoms)+header)
  countsnap=1
  coords=zeros((int(numatoms),4),dtype=float)
  sys.stdout.write("Writing [%d]: "%(snaps))
  sys.stdout.flush()
  read_atoms=1
  for line in lines:
    if line.find('ITEM: TIMESTEP')==0:
      read_atom_flag=False
      sys.stdout.write("%d "%(countsnap))
      sys.stdout.flush()
      fout.writelines("%s\nAtoms\n"%(numatoms))
      countsnap+=1
      continue
    if line.find('ITEM: ATOMS')==0:
      read_atom_flag=True
      continue
    if read_atom_flag==True:
      read_atoms+=1
      parse=line.split()
      if parse[0]!="":
        coords[int(parse[0])-1][0]=int(parse[1])
        coords[int(parse[0])-1][1]=float(parse[2])
        coords[int(parse[0])-1][2]=float(parse[3])
        coords[int(parse[0])-1][3]=float(parse[4])
      if read_atoms==int(numatoms):
        read_atoms=0
        for i in range(int(numatoms)):
          fout.writelines("%d %2.4f %2.4f %2.4f\n"%(coords[i][0],coords[i][1],coords[i][2],coords[i][3]))

  print "\nDone converting to xyz!!\n"
  fin.close()
  fout.close()
  return

if __name__ == '__main__':

    # if no input, print help and exit
    if len(sys.argv) < 2:
        print Info()
        sys.exit(1)

    inputfile=sys.argv[1]
    outfile=sys.argv[2]

    lmp2xyz(inputfile,outfile.split()[0])

