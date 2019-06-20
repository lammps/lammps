#!/usr/bin/env python

from __future__ import division

import numpy as np
import os

def read_volumes(fn):
    '''Read the volumes from a LAMMPS x y z dump'''
    def bounds2volume(bounds):
        '''Convert bounding box to cell box, then compute volume from it'''
        xy,xz,yz = bounds[0,2],bounds[1,2],bounds[2,2]
        xlo = bounds[0,0] - np.amin([0.0,xy,xz,xy+xz])
        xhi = bounds[0,1] - np.amax([0.0,xy,xz,xy+xz])
        ylo = bounds[1,0] - np.amin([0.0,yz])
        yhi = bounds[1,1] - np.amax([0.0,yz])
        zlo = bounds[2,0]
        zhi = bounds[2,1]
        # Cell box is upper triangular, so volume is product of diagonal elements
        return (xhi-xlo)*(yhi-ylo)*(zhi-zlo)
    volumes = []
    with open(fn,'r') as f:
        while True:
            line = f.readline()
            if not line: break
            if line.startswith("ITEM: TIMESTEP"):
                timestep = int(f.readline())
            if line.startswith("ITEM: BOX BOUNDS"):
                bounds = np.array([[float(w) for w in f.readline().split()] for i in range(3)])
                volume = bounds2volume(bounds)
                volumes.append([timestep, volume])
    return np.asarray(volumes)

if __name__=='__main__':
    # The following file should contain the positions obtained from an NPT run
    npt_dir = '../mil53al'
    fn_positions = os.path.join(npt_dir, 'positions.dat')
    volumes = read_volumes(fn_positions)
    # The following file is a template input for the fixed volume simulations
    # Input values that need to be adapted should be in the format __value__
    fn_template = 'template.in'
    with open(fn_template,'r') as f:
        template = f.read()
    # Request simulations at the following fixed volumes
    fixed_volumes = np.arange(1500.0,3200.0,step=25.0)
    for iv, v in enumerate(fixed_volumes):
        index = np.argmin(np.abs(volumes[:,1]-v))
        scale = np.power(volumes[index,1]/v,-1.0/3.0)
        # Check that there is a snapshot in the NPT simulation that is
        # more or less near to the requested volume
        assert np.abs(scale-1.0)<5e-2, "V wanted = %6.2f V found = %6.2f  Scale = %f" % (v/angstrom**3,allvolumes[index]/angstrom**3,scale)
        dn = 'nvst_%05d'%(iv)
        if not os.path.isdir(dn): os.makedirs(dn)
        with open(os.path.join(dn,'lammps.in'),'w') as f:
            lines = template.replace('__scale__','%8.6f'%scale)
            lines = lines.replace('__nptdir__',os.path.join('..',npt_dir))
            lines = lines.replace('__step__','%d'%(volumes[index,0]))
            f.write(lines)
