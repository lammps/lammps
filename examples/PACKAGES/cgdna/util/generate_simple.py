# Setup tool for oxDNA input in LAMMPS format.

# for python2/3 compatibility
from __future__ import print_function

import math,numpy as np,sys,os

# system size
lxmin = -115.0
lxmax = +115.0
lymin = -115.0
lymax = +115.0
lzmin = -115.0
lzmax = +115.0

# rise in z-direction
r0 = 0.7

# definition of single untwisted strand
def single():

  strand = inp[1].split(':')

  com_start=strand[0].split(',')

  posx=float(com_start[0])
  posy=float(com_start[1])
  posz=float(com_start[2])
  risex=0
  risey=0
  risez=r0

  strandstart=len(nucleotide)+1

  for letter in strand[1]:
    temp=[]

    temp.append(nt2num[letter])
    temp.append([posx,posy,posz])
    vel=[0,0,0,0,0,0]
    temp.append(vel)
    temp.append(shape)

    quat=[1,0,0,0]
    temp.append(quat)

    posx=posx+risex
    posy=posy+risey
    posz=posz+risez

    if (len(nucleotide)+1 > strandstart):
      topology.append([1,len(nucleotide),len(nucleotide)+1])

    nucleotide.append(temp)

  return

# definition of single twisted strand
def single_helix():

  strand = inp[1].split(':')

  com_start=strand[0].split(',')
  twist=0.6

  posx = float(com_start[0])
  posy = float(com_start[1])
  posz = float(com_start[2])
  risex=0
  risey=0
  risez=math.sqrt(r0**2-4.0*math.sin(0.5*twist)**2) 

  dcomh=0.76
  axisx=dcomh + posx
  axisy=posy

  strandstart=len(nucleotide)+1
  quat=[1,0,0,0]

  qrot0=math.cos(0.5*twist)
  qrot1=0
  qrot2=0
  qrot3=math.sin(0.5*twist)

  for letter in strand[1]:
    temp=[]

    temp.append(nt2num[letter])
    temp.append([posx,posy,posz])
    vel=[0,0,0,0,0,0]
    temp.append(vel)
    temp.append(shape)

    temp.append(quat)

    quat0 = quat[0]*qrot0 - quat[1]*qrot1 - quat[2]*qrot2 - quat[3]*qrot3 
    quat1 = quat[0]*qrot1 + quat[1]*qrot0 + quat[2]*qrot3 - quat[3]*qrot2 
    quat2 = quat[0]*qrot2 + quat[2]*qrot0 + quat[3]*qrot1 - quat[1]*qrot3 
    quat3 = quat[0]*qrot3 + quat[3]*qrot0 + quat[1]*qrot2 + quat[2]*qrot1 

    quat = [quat0,quat1,quat2,quat3]

    posx=axisx - dcomh*(quat[0]**2+quat[1]**2-quat[2]**2-quat[3]**2)
    posy=axisy - dcomh*(2*(quat[1]*quat[2]+quat[0]*quat[3]))
    posz=posz+risez

    if (len(nucleotide)+1 > strandstart):
      topology.append([1,len(nucleotide),len(nucleotide)+1])

    nucleotide.append(temp)

  return

# definition of twisted duplex  
def duplex():

  strand = inp[1].split(':')

  com_start=strand[0].split(',')
  twist=0.6

  compstrand=[]
  comptopo=[]

  posx1 = float(com_start[0])
  posy1 = float(com_start[1])
  posz1 = float(com_start[2])

  risex=0
  risey=0
  risez=math.sqrt(r0**2-4.0*math.sin(0.5*twist)**2) 

  dcomh=0.76
  axisx=dcomh + posx1
  axisy=posy1

  posx2 = axisx + dcomh  
  posy2 = posy1
  posz2 = posz1

  strandstart=len(nucleotide)+1

  quat1=[1,0,0,0]
  quat2=[0,0,-1,0]

  qrot0=math.cos(0.5*twist)
  qrot1=0
  qrot2=0
  qrot3=math.sin(0.5*twist)

  for letter in strand[1]:
    temp1=[]
    temp2=[]

    temp1.append(nt2num[letter])
    temp2.append(compnt2num[letter])

    temp1.append([posx1,posy1,posz1])
    temp2.append([posx2,posy2,posz2])

    vel=[0,0,0,0,0,0]
    temp1.append(vel)
    temp2.append(vel)

    temp1.append(shape)
    temp2.append(shape)

    temp1.append(quat1)
    temp2.append(quat2)

    quat1_0 = quat1[0]*qrot0 - quat1[1]*qrot1 - quat1[2]*qrot2 - quat1[3]*qrot3 
    quat1_1 = quat1[0]*qrot1 + quat1[1]*qrot0 + quat1[2]*qrot3 - quat1[3]*qrot2 
    quat1_2 = quat1[0]*qrot2 + quat1[2]*qrot0 + quat1[3]*qrot1 - quat1[1]*qrot3 
    quat1_3 = quat1[0]*qrot3 + quat1[3]*qrot0 + quat1[1]*qrot2 + quat1[2]*qrot1 

    quat1 = [quat1_0,quat1_1,quat1_2,quat1_3]

    posx1=axisx - dcomh*(quat1[0]**2+quat1[1]**2-quat1[2]**2-quat1[3]**2)
    posy1=axisy - dcomh*(2*(quat1[1]*quat1[2]+quat1[0]*quat1[3]))
    posz1=posz1+risez

    quat2_0 = quat2[0]*qrot0 - quat2[1]*qrot1 - quat2[2]*qrot2 + quat2[3]*qrot3 
    quat2_1 = quat2[0]*qrot1 + quat2[1]*qrot0 - quat2[2]*qrot3 - quat2[3]*qrot2 
    quat2_2 = quat2[0]*qrot2 + quat2[2]*qrot0 + quat2[3]*qrot1 + quat2[1]*qrot3 
    quat2_3 =-quat2[0]*qrot3 + quat2[3]*qrot0 + quat2[1]*qrot2 + quat2[2]*qrot1 

    quat2 = [quat2_0,quat2_1,quat2_2,quat2_3]

    posx2=axisx + dcomh*(quat1[0]**2+quat1[1]**2-quat1[2]**2-quat1[3]**2)
    posy2=axisy + dcomh*(2*(quat1[1]*quat1[2]+quat1[0]*quat1[3]))
    posz2=posz1

    if (len(nucleotide)+1 > strandstart):
      topology.append([1,len(nucleotide),len(nucleotide)+1])
      comptopo.append([1,len(nucleotide)+len(strand[1]),len(nucleotide)+len(strand[1])+1])

    nucleotide.append(temp1)
    compstrand.append(temp2)

  for ib in range(len(compstrand)):
    nucleotide.append(compstrand[len(compstrand)-1-ib])

  for ib in range(len(comptopo)):
    topology.append(comptopo[ib])

  return

# definition of array of duplexes  
def duplex_array():

  strand = inp[1].split(':')
  number=strand[0].split(',')
  posz1_0 = float(strand[1])
  twist=0.6

  nx = int(number[0])
  ny = int(number[1])

  dx = (lxmax-lxmin)/nx
  dy = (lymax-lymin)/ny

  risex=0
  risey=0
  risez=math.sqrt(r0**2-4.0*math.sin(0.5*twist)**2) 
  dcomh=0.76

  for ix in range(nx):

    axisx=lxmin + dx/2 + ix * dx

    for iy in range(ny):

      axisy=lymin + dy/2 + iy * dy

      compstrand=[]
      comptopo=[]

      posx1 = axisx - dcomh
      posy1 = axisy
      posz1 = posz1_0

      posx2 = axisx + dcomh  
      posy2 = posy1
      posz2 = posz1

      strandstart=len(nucleotide)+1
      quat1=[1,0,0,0]
      quat2=[0,0,-1,0]

      qrot0=math.cos(0.5*twist)
      qrot1=0
      qrot2=0
      qrot3=math.sin(0.5*twist)

      for letter in strand[2]:
        temp1=[]
        temp2=[]

        temp1.append(nt2num[letter])
        temp2.append(compnt2num[letter])

        temp1.append([posx1,posy1,posz1])
        temp2.append([posx2,posy2,posz2])

        vel=[0,0,0,0,0,0]
        temp1.append(vel)
        temp2.append(vel)

        temp1.append(shape)
        temp2.append(shape)

        temp1.append(quat1)
        temp2.append(quat2)

        quat1_0 = quat1[0]*qrot0 - quat1[1]*qrot1 - quat1[2]*qrot2 - quat1[3]*qrot3
        quat1_1 = quat1[0]*qrot1 + quat1[1]*qrot0 + quat1[2]*qrot3 - quat1[3]*qrot2
        quat1_2 = quat1[0]*qrot2 + quat1[2]*qrot0 + quat1[3]*qrot1 - quat1[1]*qrot3
        quat1_3 = quat1[0]*qrot3 + quat1[3]*qrot0 + quat1[1]*qrot2 + quat1[2]*qrot1

        quat1 = [quat1_0,quat1_1,quat1_2,quat1_3]

        posx1=axisx - dcomh*(quat1[0]**2+quat1[1]**2-quat1[2]**2-quat1[3]**2)
        posy1=axisy - dcomh*(2*(quat1[1]*quat1[2]+quat1[0]*quat1[3]))
        posz1=posz1+risez

        quat2_0 = quat2[0]*qrot0 - quat2[1]*qrot1 - quat2[2]*qrot2 + quat2[3]*qrot3
        quat2_1 = quat2[0]*qrot1 + quat2[1]*qrot0 - quat2[2]*qrot3 - quat2[3]*qrot2
        quat2_2 = quat2[0]*qrot2 + quat2[2]*qrot0 + quat2[3]*qrot1 + quat2[1]*qrot3
        quat2_3 =-quat2[0]*qrot3 + quat2[3]*qrot0 + quat2[1]*qrot2 + quat2[2]*qrot1

        quat2 = [quat2_0,quat2_1,quat2_2,quat2_3]

        posx2=axisx + dcomh*(quat1[0]**2+quat1[1]**2-quat1[2]**2-quat1[3]**2)
        posy2=axisy + dcomh*(2*(quat1[1]*quat1[2]+quat1[0]*quat1[3]))
        posz2=posz1

        if (len(nucleotide)+1 > strandstart):
          topology.append([1,len(nucleotide),len(nucleotide)+1])
          comptopo.append([1,len(nucleotide)+len(strand[2]),len(nucleotide)+len(strand[2])+1])

        nucleotide.append(temp1)
        compstrand.append(temp2)

      for ib in range(len(compstrand)):
        nucleotide.append(compstrand[len(compstrand)-1-ib])

      for ib in range(len(comptopo)):
        topology.append(comptopo[ib])

  return

# main part
nt2num = {'A':1, 'C':2, 'G':3, 'T':4}
compnt2num = {'T':1, 'G':2, 'C':3, 'A':4}
shape = [1.1739845031423408,1.1739845031423408,1.1739845031423408]

nucleotide=[]
topology=[]

seqfile = open(sys.argv[1],'r')

# process sequence file line by line
for line in seqfile:

  inp = line.split()
  if inp[0] == 'single':
    single()
  if inp[0] == 'single_helix':
    single_helix()
  if inp[0] == 'duplex':
    duplex()
  if inp[0] == 'duplex_array':
    duplex_array()

# output atom data in LAMMPS format
out = open(sys.argv[2],'w')

out.write('# LAMMPS data file\n')
out.write('%d atoms\n' % len(nucleotide))
out.write('%d ellipsoids\n' % len(nucleotide))
out.write('%d bonds\n' % len(topology))
out.write('\n')
out.write('4 atom types\n')
out.write('1 bond types\n')
out.write('\n')
out.write('# System size\n')
out.write('%f %f xlo xhi\n' % (lxmin,lxmax))
out.write('%f %f ylo yhi\n' % (lymin,lymax))
out.write('%f %f zlo zhi\n' % (lzmin,lzmax))
out.write('\n')
out.write('Masses\n')
out.write('\n')
out.write('1 3.1575\n')
out.write('2 3.1575\n')
out.write('3 3.1575\n')
out.write('4 3.1575\n')

out.write('\n')
out.write('# Atom-ID, type, position, molecule-ID, ellipsoid flag, density\n')
out.write('Atoms\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write("%d %d %22.16le %22.16le %22.16le 1 1 1\n" % (ib+1,nucleotide[ib][0],nucleotide[ib][1][0],nucleotide[ib][1][1],nucleotide[ib][1][2]))

out.write('\n')
out.write('# Atom-ID, translational, rotational velocity\n')
out.write('Velocities\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write("%d %22.16le %22.16le %22.16le %22.16le %22.16le %22.16le\n" % (ib+1,nucleotide[ib][2][0],nucleotide[ib][2][1],nucleotide[ib][2][2],nucleotide[ib][2][3],nucleotide[ib][2][4],nucleotide[ib][2][5]))

out.write('\n')
out.write('# Atom-ID, shape, quaternion\n')
out.write('Ellipsoids\n')
out.write('\n')
for ib in range(len(nucleotide)):
  out.write("%d %22.16le %22.16le %22.16le %22.16le %22.16le %22.16le %22.16le\n" % (ib+1,nucleotide[ib][3][0],nucleotide[ib][3][1],nucleotide[ib][3][2],nucleotide[ib][4][0],nucleotide[ib][4][1],nucleotide[ib][4][2],nucleotide[ib][4][3]))

out.write('\n')
out.write('# Bond topology\n')
out.write('Bonds\n')
out.write('\n')
for ib in range(len(topology)):
  out.write("%d %d %d %d\n" % (ib+1,topology[ib][0],topology[ib][1],topology[ib][2]))

out.close() 

seqfile.close()
sys.exit(0)


