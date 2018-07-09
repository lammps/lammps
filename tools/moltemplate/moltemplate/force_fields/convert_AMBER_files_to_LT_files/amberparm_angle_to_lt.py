#!/usr/bin/env python

import sys

lines_gaff = sys.stdin.readlines()
angle_style_name = 'harmonic'

sys.stdout.write('  write_once("In Settings") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    atypes = line[:8].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    at3 = atype3.replace('X','*')
    angletype = '@angle:'+atype1+'-'+atype2+'-'+atype3

    tokens= line[8:].split()
    keq = tokens[0]
    req = tokens[1]
    comments=' '.join(tokens[2:])
    sys.stdout.write('    angle_coeff '+angletype+' '+angle_style_name+' '+keq+' '+req)
    if len(comments.strip()) > 0:
        sys.stdout.write('   # '+comments)
    sys.stdout.write('\n')


sys.stdout.write('  } # (end of angle_coeffs)\n')
sys.stdout.write('\n')
sys.stdout.write('  write_once("Data Angles By Type") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    atypes = line[:8].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    atype3 = atypes[2].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    at3 = atype3.replace('X','*')
    angletype = '@angle:'+atype1+'-'+atype2+'-'+atype3

    #tokens= line[8:].split()
    #keq = tokens[0]
    #req = tokens[1]
    #comments=' '.join(tokens[2:])
    sys.stdout.write('    '+angletype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+'\n')

sys.stdout.write('  } # (end of Angles By Type)\n')
sys.stdout.write('\n')
