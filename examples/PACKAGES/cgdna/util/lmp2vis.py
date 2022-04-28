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

Program: lmp2vis.py

Produces a simple representation of the oxDNA nucleotide with separate 
particles for backbone and base interaction site. The base particle inherits 
the atom type, whereas the backbone particle acquires an offset of 10. This
can be changed below.

Usage: 
$$ python lmp2vis.py [visualisation (vmd OR ovito, default=ovito)] input_filename output_filename

Requirements:
The LAMMPS trajectory input file needs to contain the following data columns:
id mol type x y z vx vy vz c_quat[1] c_quat[2] c_quat[3] c_quat[4]
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
    typb = typ + 10  # defines new backbone types, here offset is 10 from atom (base) type
    x, y, z = float(list1[colind[3]]), float(list1[colind[4]]), float(list1[colind[5]])
    vx, vy, vz = float(list1[colind[6]]), float(list1[colind[7]]), float(list1[colind[8]])
    c_quat1, c_quat2, c_quat3, c_quat4 = \
            float(list1[colind[9]]), float(list1[colind[10]]), float(list1[colind[11]]), float(list1[colind[12]])

    ex, ey, ez = q_to_exyz(c_quat1, c_quat2, c_quat3, c_quat4)

    # position of sugar-phosphate backbone interaction site in oxDNA2
    x1, y1, z1 = x -0.34*ex[0]+0.3408*ey[0],y -0.34*ex[1]+0.3408*ey[1], z-0.34*ex[2]+0.3408*ey[2]

    # position of base interaction site in oxDNA/oxDNA2
    x2, y2, z2 = x +0.4*ex[0], y + 0.4*ex[1], z+0.4*ex[2]

    # compose basic output data: id, molecule id, type, position, velocity quaternion
    line1 = '%d'%(2*ident-1) +' '+ '%d'%mol +' '+ '%d'%typb +' '+\
        '%13.6e'%(x1) +' '+ '%13.6e'%(y1) +' '+ '%13.6e'%(z1) +' '+\
        '%13.6e'%(vx) +' '+ '%13.6e'%(vy) +' '+ '%13.6e'%(vz) +' '+\
        '%13.6e'%(c_quat1) +' '+ '%13.6e'%(c_quat2) +' '+ '%13.6e'%(c_quat3) +' '+ '%13.6e'%(c_quat4)

    line2 = '%d'%(2*ident) +' '+ '%d'%mol +' '+ '%d'%typ +' '+\
        '%13.6e'%(x2) +' '+ '%13.6e'%(y2) +' '+ '%13.6e'%(z2) +' '+\
        '%13.6e'%(vx) +' '+ '%13.6e'%(vy) +' '+ '%13.6e'%(vz) +' '+\
        '%13.6e'%(c_quat1) +' '+ '%13.6e'%(c_quat2) +' '+' %13.6e'%(c_quat3) +' '+ '%13.6e'%(c_quat4)

    # use oblate particles for bases in ovito
    shape_sphere = ' 0.4 0.4 0.4'
    shape_ellipsoid = ' 0.5 0.2 0.1'
    if vismethod == 'ovito':
        line1 += shape_sphere +' '
        line2 += shape_ellipsoid +' '

    line1 += '\n'
    line2 += '\n'

    line = line1 + line2
    return line

### main part ###

# digest command line input
if len(sys.argv)<3:
    print("Syntax: $$ python lmp2vis.py [visualisation (vmd OR ovito, default=ovito)] input_filename output_filename")
    sys.exit(1)

if len(sys.argv)==3:
    vismethod = 'ovito' # default visualisation method
    infilename  = sys.argv[1]
    outfilename = sys.argv[2]

if len(sys.argv)==4:
    vismethod = sys.argv[1]
    if (sys.argv[1]!='vmd' and sys.argv[1]!='ovito'):
        vismethod = 'ovito' # default visualisation method
    infilename  = sys.argv[2]
    outfilename = sys.argv[3]

print('# Converting LAMMPS output for visualisation with %s' % vismethod)     

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

    # find number of atoms in timestep
    if line.find('ITEM: NUMBER OF ATOMS') != -1: 
        w.write(line)
        N=int(r.readline())
        # write to output file and read next line
        w.write('%d'%int(2*N)+'\n')
        line=r.readline()
        n+=2

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
            vxindex = linestring.index('vx')
            vyindex = linestring.index('vy')
            vzindex = linestring.index('vz')
            qwindex = linestring.index('c_quat[1]')
            qxindex = linestring.index('c_quat[2]')
            qyindex = linestring.index('c_quat[3]')
            qzindex = linestring.index('c_quat[4]')

            # create header
            header = linestring[0] + ' ' + linestring[1] + ' ' + \
                    linestring[idindex] + ' ' + linestring[molindex]+ ' ' + linestring[typeindex]+ ' ' + \
                    linestring[xindex]+ ' ' + linestring[yindex]+ ' ' + linestring[zindex]+ ' ' + \
                    linestring[vxindex]+ ' ' + linestring[vyindex]+ ' ' + linestring[vzindex]+ ' ' + \
                    linestring[qwindex]+ ' ' + linestring[qxindex]+ ' ' + linestring[qyindex]+ ' ' + linestring[qzindex]

            # extend header for ovito
            if vismethod == 'ovito':
                header += ' shape[0] shape[1] shape[2]'
            header += '\n'

            # store column number in data line, -2 offset form header line
            colind = [idindex-2,molindex-2,typeindex-2,xindex-2,yindex-2,zindex-2,\
                      vxindex-2,vyindex-2,vzindex-2,qwindex-2,qxindex-2,qyindex-2,qzindex-2]
            pass1 = 1

        # begin processing atom data in timestep
        w.write(header)
        # tranform each atom and write to output file
        for i in range(N):
            line=r.readline()
            w.write(transform(line,colind))
            n+=1
        # end processing atom data

    # duplicate all lines that are not in any other section 
    else:
        w.write(line)
        n+=1

print('# Done                                ')     

r.close()
w.close()
