#!/usr/bin/env python
# check for inconsistent thermodata (energy/pressure)
# second line of the first run has to match first line
# each of the second and third. second line of second
# and third run have to match as well.
# Example output:
#############################################################################
#Step Temp E_pair E_mol TotEng Press 
#       0            0            0  0.040304848  0.040304848      183.537503 
#       1   0.22393239            0  0.036893575  0.040231075      175.827611 
#--
#Step Temp E_pair E_mol TotEng Press 
#       0   0.22393239            0  0.036893575  0.040231075      175.827611 
#       1   0.81991975            0  0.027814603  0.040034728      153.287062 
#--
#Step Temp E_pair E_mol TotEng Press 
#       1   0.22393239            0  0.036893575  0.040231075      175.827611 
#       2   0.81991975            0  0.027814603  0.040034728      153.287062 
#############################################################################

import sys

for j in range(1,len(sys.argv)):
    f = sys.argv[j]
    fp = open(f,'r')
    fp.readline()
    fp.readline()
    line = fp.readline()
    col1 = line.split()
    fp.readline()
    fp.readline()
    line = fp.readline()
    col2 = line.split()
    line = fp.readline()
    col3 = line.split()
    fp.readline()
    fp.readline()
    line = fp.readline()
    col4 = line.split()
    line = fp.readline()
    col5 = line.split()
    fp.close()
    diff = 0
    for i in [1,3,4,5]:
        if (col1[i] != col2[i]):
            diff = 1
        if (col1[i] != col4[i]):
            diff = 1
        if (col3[i] != col5[i]):
            diff = 1
    if (diff > 0):
        print "Inconsistent thermo data"
        exit(1)

