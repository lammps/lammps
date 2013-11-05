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
    impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

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
    impropertype = '@improper:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

    sys.stdout.write('    '+impropertype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')

sys.stdout.write('  } # (end of Impropers By Type)\n')
sys.stdout.write('\n')

