#!/usr/bin/env python

import sys

lines_gaff = sys.stdin.readlines()

#pair_style = 'lj/charmm/coul/long'

    # NOTE: Long-range coulombic forces were disabled intentionally. (See below)
    #       If you want to use long-range electrostatics, uncomment these lines:
    # Instead I use hybrid lj/charmm/coul/charmm by default, because
    # LAMMPS complains if you attempt to use lj/charmm/coul/long on a
    # system if it does not contain any charged particles.
    # Currently, moltemplate does not assign atomic charge, 
    # so this problem occurs frequently.

#pair_style = 'lj/charmm/coul/charmm'
pair_style = 'lj/charmm/coul/long'

sys.stdout.write('  write_once(\"In Settings\") {\n')

for i in range(0, len(lines_gaff)):
    line = lines_gaff[i]
    tokens= line.split()
    atype = tokens[0]
    sig=tokens[1]
    eps=tokens[2]
    comments=' '.join(tokens[3:])
    sys.stdout.write('    pair_coeff @atom:'+atype+' @atom:'+atype+' '+pair_style+' '+eps+' '+sig+'   # '+comments+'\n')

sys.stdout.write('  } # (end of pair_coeffs)\n')
sys.stdout.write('\n')
