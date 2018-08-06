#!/usr/bin/env python

import sys

lines_gaff = sys.stdin.readlines()
bond_style_name = 'harmonic'

sys.stdout.write('  write_once("In Settings") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    tokens= line.split()
    atypes = line[:6].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    bondtype = '@bond:'+atype1+'-'+atype2

    tokens= line[5:].split()
    keq = tokens[0]
    req = tokens[1]
    comments=' '.join(tokens[2:])
    sys.stdout.write('    bond_coeff '+bondtype+' '+bond_style_name+' '+keq+' '+req)
    if len(comments.strip()) > 0:
        sys.stdout.write('   # '+comments)
    sys.stdout.write('\n')


sys.stdout.write('  } # (end of bond_coeffs)\n')
sys.stdout.write('\n')
sys.stdout.write('  write_once("Data Bonds By Type") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    atypes = line[:6].split('-')
    atype1 = atypes[0].strip()
    atype2 = atypes[1].strip()
    at1 = atype1.replace('X','*')
    at2 = atype2.replace('X','*')
    bondtype = '@bond:'+atype1+'-'+atype2

    #tokens= line[5:].split()
    #keq = tokens[0]
    #req = tokens[1]
    #comments=' '.join(tokens[2:])
    sys.stdout.write('    '+bondtype+' @atom:'+at1+' @atom:'+at2+'\n')

sys.stdout.write('  } # (end of Bonds By Type)\n')
sys.stdout.write('\n')
