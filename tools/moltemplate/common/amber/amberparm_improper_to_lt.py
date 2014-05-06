#!/usr/bin/env python

import sys

lines_gaff = sys.stdin.readlines()
improper_style_name = 'cvff'

sys.stdout.write('  write_once("In Settings") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    atypes = line[:11].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    atype4 = atypes[3].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    at3 = atype3.replace('X','*')
    at4 = atype4.replace('X','*')
    #impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4
    #sys.stdout.write('    '+impropertype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')
    # Oops.  This is incorrect.
    # In moltemplate, the central atom is the first atom,
    # In "gaff.dat", the central atom is the third atom
    # http://archive.ambermd.org/201307/0519.html
    impropertype = '@improper:'+atype3+'-'+atype1+'-'+atype2+'-'+atype4

    tokens= line[11:].split()
    Kn = float(tokens[0])
    dn = float(tokens[1])
    n = int(float(tokens[2]))
    comments=' '.join(tokens[3:])

    if (dn < 0.001):
        sys.stdout.write('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' 1 '+str(n)+'   # '+comments+'\n')
    elif (179.999 < abs(dn) < 180.001):
        sys.stdout.write('    improper_coeff '+impropertype+' '+improper_style_name+' '+str(Kn)+' -1 '+str(n)+'   # '+comments+'\n')
    else:
        sys.stderr.write('Error: Illegal bondImproper parameters:\n'
                         '       As of 2013-8-03, LAMMPS doens hot have an improper style\n'
                         '       which can handle impropers with gamma != 0 or 180\n')
        exit(-1)



sys.stdout.write('  } # (end of improper_coeffs)\n')
sys.stdout.write('\n')
sys.stdout.write('  write_once("Data Impropers By Type") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    atypes = line[:11].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    atype4 = atypes[3].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    at3 = atype3.replace('X','*')
    at4 = atype4.replace('X','*')

    #impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4
    #sys.stdout.write('    '+impropertype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')
    # Oops.  This is incorrect.
    # In moltemplate, the central atom is the first atom,
    # In "gaff.dat", the central atom is the third atom
    # http://archive.ambermd.org/201307/0519.html
    impropertype = '@improper:'+atype3+'-'+atype1+'-'+atype2+'-'+atype4
    sys.stdout.write('    '+impropertype+' @atom:'+at3+' @atom:'+at1+' @atom:'+at2+' @atom:'+at4+'\n')
    

sys.stdout.write('  } # (end of Impropers By Type)\n')
sys.stdout.write('\n')

# NOTE: AMBER documentation is not clear how the improper angle is defined.
#       It's not clear if we should be using the dihedral angle between
#       planes I-J-K and J-K-L.  As of 2014-4, improper_style cvff does this.
#       Even if we create improper interactions with the angle defined between
#       the wrong planes, at least the minima should be the same 
#       (0 degrees or 180 degrees).
#       So I'm not too worried we are getting this detail wrong long as 
#       we generate new impropers realizing that the 3rd atom (K) is the 
#       central atom (according to AMBER conventions).
#
# http://structbio.vanderbilt.edu/archives/amber-archive/2007/0408.php
#
# Currently, we only apply improper torsional angles for atoms 
# in a planar conformations. Is it clear?
# Junmei 
