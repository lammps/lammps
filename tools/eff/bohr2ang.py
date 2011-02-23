Info="""
Module name: bohr2ang.py 

Author: (c) Andres Jaramillo-Botero
California Institute of Technology
ajaramil@caltech.edu
Project: pEFF
Version: August 2009

Usage: python bohr2ang.py
>>Name of data file (bohr): [datafile]

Results:
creates a datafile with extension .ang in real units

"""
import os

currdir=os.getcwd()
datafile=raw_input("Name of data file (bohr): ")
bohr2ang=0.529177249
bperatu2angperfs=0.512396271120794
f=open(currdir+'/'+datafile,'r')
w=open(currdir+'/'+datafile+'.ang','w')
lines=f.readlines()
atom_flag=False
vel_flag=False
for line in lines:
  if line.find("xlo") > 0:
    parse=line.split()
    w.write("%f %f xlo xhi\n"%(float(parse[0])*bohr2ang,float(parse[1])*bohr2ang))
  elif line.find("ylo") > 0:
    parse=line.split()
    w.write("%f %f ylo yhi\n"%(float(parse[0])*bohr2ang,float(parse[1])*bohr2ang))
  elif line.find("zlo") > 0:
    parse=line.split()
    w.write("%f %f zlo zhi\n"%(float(parse[0])*bohr2ang,float(parse[1])*bohr2ang))
  elif line.find("xy") >= 0:
    parse=line.split()
    w.write("%f %f %f xy xz yz\n"%(float(parse[0])*bohr2ang,float(parse[1])*bohr2ang,float(parse[2])*bohr2ang))
  elif atom_flag and line.strip():
    parse=line.split()
    id=parse[0]
    type=parse[1]
    q=parse[2]
    spin=parse[3]
    eradius=float(parse[4])*bohr2ang
    x=float(parse[5])*bohr2ang
    y=float(parse[6])*bohr2ang
    z=float(parse[7])*bohr2ang
    rest=" ".join(parse[8:])
    w.write("%s %s %s %s %f %f %f %f %s\n"%(id,type,q,spin,eradius,x,y,z,rest))
  elif line.find("Atoms") >= 0:
    w.write(line)
    atom_flag=True
    continue
  elif vel_flag and line != "\n":
    parse=line.split()
    id=parse[0]
    vx=float(parse[1])*bperatu2angperfs
    vy=float(parse[2])*bperatu2angperfs
    vz=float(parse[3])*bperatu2angperfs
    erv=float(parse[4])*bperatu2angperfs
    w.write("%s %f %f %f\n"%(id,vx,vy,vz,erv))
  elif line.find("Velocities") >= 0:
    w.write(line)
    atom_flag=False
    vel_flag=True
    continue
  else:
    w.write(line)

f.close()
w.close()
