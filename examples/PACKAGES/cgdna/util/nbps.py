#!/usr/bin/env python
"""
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

Program:

Usage: 
$$ python 

Requirements:
The LAMMPS trajectory input file needs to contain the following data columns:
id mol type x y z c_quat[1] c_quat[2] c_quat[3] c_quat[4]
"""

import sys, math, subprocess

# converts quaternion DOF into local body reference frame
def q_to_exyz(q1,q2,q3,q4):

    q = [q1, q2, q3, q4]
    ex = [0, 0, 0]
    ey = [0, 0, 0]
    ez = [0, 0, 0]

    ex[0]=q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3]
    ex[1]=2*(q[1]*q[2]+q[0]*q[3])
    ex[2]=2*(q[1]*q[3]-q[0]*q[2])

    ey[0]=2*(q[1]*q[2]-q[0]*q[3])
    ey[1]=q[0]*q[0]-q[1]*q[1]+q[2]*q[2]-q[3]*q[3]
    ey[2]=2*(q[2]*q[3]+q[0]*q[1])

    ez[0]=2*(q[1]*q[3]+q[0]*q[2])
    ez[1]=2*(q[2]*q[3]-q[0]*q[1])
    ez[2]=q[0]*q[0]-q[1]*q[1]-q[2]*q[2]+q[3]*q[3]

    return ex,ey,ez

# processes line by line of LAMMPS trajectory
def transform(line, colind): 

    list1 = line.split()
    ident, mol, typ = int(list1[colind[0]]), int(list1[colind[1]]), int(list1[colind[2]])
    x, y, z = float(list1[colind[3]]), float(list1[colind[4]]), float(list1[colind[5]])
    c_quat1, c_quat2, c_quat3, c_quat4 = \
            float(list1[colind[6]]), float(list1[colind[7]]), float(list1[colind[8]]), float(list1[colind[9]])

    ex, ey, ez = q_to_exyz(c_quat1, c_quat2, c_quat3, c_quat4)

    # position of base interaction site in oxDNA/oxDNA2
    x2, y2, z2 = x +0.4*ex[0], y + 0.4*ex[1], z+0.4*ex[2]

    # compose basic output data: id, molecule id, type, position
    line2 = [ident, mol, typ, x2, y2, z2]

    return line2

### main part ###

# digest command line input
if len(sys.argv)<3:
    print("Syntax: $$ python nbps.py input_filename output_filename")
    sys.exit(1)

if len(sys.argv)==3:
    infilename  = sys.argv[1]
    outfilename = sys.argv[2]

print('# Calculating number of base paris from LAMMPS trajectory output')     

# count lines to process for progress report
n = 0
try:
    result = subprocess.run(['wc', '-l', '%s'%infilename], stdout=subprocess.PIPE)
    reslist=str(result).split()
    nlines=float(reslist[5])
except:
    nlines = 100

r=open(infilename,'r')
w=open(outfilename,'w+')

pass1 = 0

for line in r:

    sys.stdout.write('# Processed %3d %%\r' % (100*n/nlines))     

   # find timestep number
    if line.find('ITEM: TIMESTEP') != -1:
        t=int(r.readline())
        n+=1

   # find number of atoms in timestep
    if line.find('ITEM: NUMBER OF ATOMS') != -1:
        N=int(r.readline())
        n+=1

    # find beginning of atom data section
    if line.find('ITEM: ATOMS') != -1:
        # first pass: extract column number of ID, molecule ID, type, postion, velocity, quaternion in header line
        if pass1 == 0:
            linestring=line.split()
            idindex = linestring.index('id')
            molindex = linestring.index('mol')
            typeindex = linestring.index('type')
            xindex = linestring.index('x')
            yindex = linestring.index('y')
            zindex = linestring.index('z')
            qwindex = linestring.index('c_quat[1]')
            qxindex = linestring.index('c_quat[2]')
            qyindex = linestring.index('c_quat[3]')
            qzindex = linestring.index('c_quat[4]')

            # store column number in data line, -2 offset form header line
            colind = [idindex-2,molindex-2,typeindex-2,xindex-2,yindex-2,zindex-2,\
                      qwindex-2,qxindex-2,qyindex-2,qzindex-2]
            pass1 = 1

        # begin processing atom data
        atom = []
        for i in range(N):
            line=r.readline()
            atom.append(transform(line,colind))
            n +=1

        rc2 = 0.2
        nbp = 0
        for i in range(N): 
            for j in range(i+1,N): 
                if (atom[i][2]+atom[j][2])%4 and atom[i][1]!=atom[j][1]:
                    r2 = (atom[i][3]-atom[j][3])*(atom[i][3]-atom[j][3]) +\
                          (atom[i][4]-atom[j][4])*(atom[i][4]-atom[j][4]) +\
                          (atom[i][5]-atom[j][5])*(atom[i][5]-atom[j][5])
                    if r2<rc2:

                        nbp += 1
        w.write('%d %d\n' % (t,nbp))
        # end processing atom data

    n+=1

print('# Done                                ')     

r.close()
w.close()
