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
from math import log10,floor
from numpy import zeros

masses={"1.00794":"H","4.002602":"He","6.941":"Li","9.012182":"Be","10.811":"B","12.0107":"C","1.00":"Au","0.0005486":"Au"}
mass_floor={1:"H",4:"He",6:"Li",9:"Be",10:"B",12:"C",0:"Au",28:"Si"}

def lmp2xyz(lammps,xyz,xpos):
  print "\nGenerating %s file"%(xyz)
  fin=open(lammps,'r')
  fout=open(xyz,'w')
  data=raw_input("Do you have a corresponding data file? please enter filename or 'n': ")
  count=1
  if data!='n': 
    dataf=open(data,'r')
    datafile=dataf.readlines()
    dataf.close()
    for line in datafile:
      if line.find("atom types")>=0:
        numtypes=int(line.split()[0])
        mass=zeros(numtypes,dtype=float)
      elif line.find("Masses")>=0:
        count+=1+datafile.index(line)
      elif line.find("Atoms")>=0:
        break
    for i in range(numtypes):
      mass[i]=float(datafile[count].split()[1])
      count+=1
  else: 
    print "\nWill continue without a data file specification"
  header=9
  lines=fin.readlines()
  numatoms=lines[3].split()[0]
  fsize=os.system("wc -l %s> lines"%(lammps))
  tmp=open('lines','r')
  tlines=tmp.readline()
  tmp.close()
  os.system("rm lines")
  flines=int(tlines.split()[0])
  snaps=flines/(int(numatoms)+header)
  countsnap=1
  if data!='n': coords={}
  else: coords=zeros((int(numatoms),4),dtype=float)
#  sys.stdout.write("Writing %d snapshots\n"%(snaps))
#  sys.stdout.flush()
  read_atoms=0
  types={}
  for line in lines:
    if line.find('ITEM: TIMESTEP')==0:
      read_atom_flag=False
#      sys.stdout.write("%d "%(countsnap))
#      sys.stdout.flush()
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
        if data!='n': 
          if parse[1] not in types.keys():
            type=raw_input("Atom name for type %s: "%parse[1])
            types[parse[1]]=type
          coords[int(parse[0])-1]=[types[parse[1]],float(parse[xpos-1]),float(parse[xpos]),float(parse[xpos+1])]
        else: 
          coords[int(parse[0])-1][0]=int(parse[1])
          coords[int(parse[0])-1][1]=float(parse[xpos-1])
          coords[int(parse[0])-1][2]=float(parse[xpos])
          coords[int(parse[0])-1][3]=float(parse[xpos+1])
      if read_atoms==int(numatoms):
        read_atoms=0
        for i in range(int(numatoms)):
          if data!='n': fout.writelines("%s %2.4f %2.4f %2.4f\n"%(coords[i][0],coords[i][1],coords[i][2],coords[i][3]))
          else: fout.writelines("%d %2.4f %2.4f %2.4f\n"%(coords[i][0],coords[i][1],coords[i][2],coords[i][3]))

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
    if len(sys.argv)==4:
      xpos=sys.arv[3]-1
    else: xpos=5
    lmp2xyz(inputfile,outfile.split()[0],xpos)

