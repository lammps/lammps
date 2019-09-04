#!/usr/bin/env python

import os, sys, numpy as np

tempfn = os.path.abspath(sys.argv[1])
logfnprefix = os.path.abspath(sys.argv[2])

temps = np.loadtxt(tempfn)
ntemps = len(temps)
u_kn = []
start_token = 'Step Temp PotEng'
end_token = 'Loop time'

for i in range(ntemps):
    logfn = '%s.%d' % (logfnprefix, i)
    with open(logfn, 'r') as of:
        lines = of.readlines()
    
    # extract relevant lines from logfile
    start = [lines.index(line) for line in lines if line.startswith(start_token)][0]
    lines = lines[(start+1) : ]
    stop = [lines.index(line) for line in lines if line.startswith(end_token)][0]
    lines = lines[:stop]
    
    # store the potential energies
    pe = [float(line.strip().split()[-1]) for line in lines]
    u_kn.append(pe)
    
u_kn = np.array(u_kn)
np.savetxt('ene.peptide', u_kn, fmt = '%5.5f')
