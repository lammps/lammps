#!/usr/bin/env python

from __future__ import division

import numpy as np
import os
from glob import glob
import subprocess
from scipy import interpolate
import matplotlib.pyplot as plt
import matplotlib

def read_lammps_volumes_pressures(fn):
    '''Read time, volume and pressure from LAMMPS log file'''
    out,err = subprocess.Popen(["grep","Time",fn],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    time = []
    for line in out.decode().split('\n'):
        w = line.split()
        if len(w)!=9: continue
        time.append(float(w[5]))
    time = np.asarray(time)
    out,err = subprocess.Popen(["grep","Volume",fn],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    volumes = []
    for line in out.decode().split('\n'):
        w = line.split()
        if len(w)!=9: continue
        volumes.append(float(w[-1]))
    volumes = np.asarray(volumes)
    out,err = subprocess.Popen(["grep","Press",fn],stdout=subprocess.PIPE,stderr=subprocess.PIPE).communicate()
    pressures = []
    for line in out.decode().split('\n'):
        w = line.split()
        if len(w)!=3: continue
        pressures.append(float(w[2]))
    pressures = np.asarray(pressures)
    imax = np.amin([pressures.shape[0], volumes.shape[0], time.shape[0]])
    return time[:imax], volumes[:imax], pressures[:imax]

def thermo_integration_spline(v0,p0):
    '''
    Integrate the pressure over the volume to obtain the free energy.
    '''
    # Set the spacing of the interior knots:
    # low values result in a spline that follows data closely, but might be noisy
    # high values result in a spline that is smooth, but might deviate from data
    window = 5
    # Construct the spline
    spline = interpolate.LSQUnivariateSpline(v0, p0, v0[1:-1:window])
    # Find roots with negative slope
    dp_spline = spline.derivative()
    roots = spline.roots()
    roots = roots[dp_spline(roots)<0.0]
#    dp_spline2 = interpolate.LSQUnivariateSpline(v0, dp_spline(v0), v0[1:-1:window])
    assert roots.shape[0] == 1 
    trans_pressures = np.amax(p0)
    # Evaluate spline at volume points
    p0_spline = spline(v0)
    # Integrate spline
    spline_int = spline.antiderivative()
    # Evaluate minus the integral at volume points, this is the free energy
    f0_spline = -spline_int(v0)
    # Shift free energy so minimum equals zero
    f0_spline -= np.amin(f0_spline)
    f_minima = -spline_int(roots)
    f_minima -= np.amin(f_minima)
    return f0_spline, p0_spline, roots, f_minima, trans_pressures

if __name__=='__main__':
    # The same units as in the LAMMPS scripts are retained, in this case real
    # units. For example time in femtoseconds, volume in A**3, pressure in atm, ...
    teq = 2500 # Only include samples beyond teq for calculating averages
    # Collect data from all fixed-volume simulations
    pattern = 'nvst*/lammps.out'
    fns = sorted(glob(pattern))
    data = []
    for fn in fns:
        dn = os.path.dirname(fn)
        time, volume, pressure = read_lammps_volumes_pressures(fn)
        mask = time>teq
        if np.sum(mask)==0: continue
        tmax = time[-1]
        N = np.sum(mask)
        V = np.mean(volume[mask])
        # Check that the volume was constant
        assert np.std(volume[mask]) < 1e-4, "%40s %5d %8.2e"%(fn,N,np.std(volume[mask]))
        P = np.mean(pressure[mask])
        Perr = np.std(pressure[mask])
        row = [N,tmax,V,P,Perr]
        print("%12s Nsamples = %6d tmax = %8.1f fs V = %8.2f A^3 P = %8.2f (+-%8.2f) atm" % (fn,row[0],row[1],row[2],row[3],row[4]))
        data.append([row[2],row[3]])
    data = np.asarray(data)
    # Perform thermodynamic integration to get the free energy as a function
    # of the volume
    f, pf, roots, fmin, p_tr = thermo_integration_spline(data[:,0], data[:,1])
    plt.clf()
    plt.subplot(2,1,1)
    plt.plot(data[:,0],data[:,1],marker='o')
    plt.xlabel("Volume [A**3]")
    plt.ylabel("Pressure [atm]")
    plt.subplot(2,1,2)
    # Conversion from atm*A**3 to kJ/mol
    conversion = 6.101934874875e-05
    plt.plot(data[:,0],f*conversion,marker='o')
    plt.xlabel("Volume [A**3]")
    plt.ylabel("Free energy [kJ/mol]")
    plt.tight_layout()
    plt.savefig('profile.png')
