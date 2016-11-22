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

    # UGGHHH 

    # OLD CODE:
    #sig=tokens[1]

    #   CORRECTION #1
    # It looks the number in this part of the file is an atom radii, not a 
    # diameter.  In other words, this number is 0.5*sigma instead of sigma.
    # So we multiply it by 2.0.
    #sig=str(2.0*float(tokens[1]))
    #
    #   CORRECTION #2
    # It also appears as though they are using this convention for LennardJones
    # U(r)=epsilon*((s/r)^12-2*(s/r)^6)    instead of  4*eps*((s/r)^12-(s/r)^6)
    #    ...where "s" is shorthand for "sigma"..
    # This means we must ALSO multiply sigma in gaff.dat by 2**(-1.0/6)
    # (This change makes the two U(r) formulas equivalent.)

    #  I had to figure this out by iterations of trial and error.
    #  The official AMBER documentation is quite vague about the LJ parameters.
    #  My apologies to everyone effected by this bug!  -Andrew 2014-5-19
    # http://ambermd.org/formats.html#parm.dat
    # http://structbio.vanderbilt.edu/archives/amber-archive/2009/5072.php)

    sig=str(float(tokens[1])*2.0*pow(2.0, (-1.0/6.0)))
    eps=tokens[2]
    comments=' '.join(tokens[3:])
    sys.stdout.write('    pair_coeff @atom:'+atype+' @atom:'+atype+' '+pair_style+' '+eps+' '+sig+'   # '+comments+'\n')

sys.stdout.write('  } # (end of pair_coeffs)\n')
sys.stdout.write('\n')
