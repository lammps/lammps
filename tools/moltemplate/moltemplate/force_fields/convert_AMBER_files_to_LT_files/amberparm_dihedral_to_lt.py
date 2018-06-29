#!/usr/bin/env python

# SOME UGLY CODE HERE

import sys

lines_gaff = sys.stdin.readlines()
dihedral_style_name = 'fourier'
in_dihedral_coeffs = []


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
    dihedraltype = '@dihedral:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

    tokens= line[11:].split()
    npth = float(tokens[0])
    Kn = float(tokens[1])
    Kn /= npth  # The coeff for each fourier term is Kn/npth
                # ...I THINK (?).  (Very confusing.  See documentation below...)
    dn = float(tokens[2])
    n = int(float(tokens[3]))
    comments=' '.join(tokens[4:])
    if len(comments.strip()) > 0:
        comments = '    # ' + comments
    in_dihedral_coeffs.append([dihedraltype, Kn, n, dn, comments])
    #print(Kn, n, dn)

#for entry in in_dihedral_coeffs:
#    print(entry)
#exit()


# ---- processing dihedral fourier series ----
# ---- (negative "n" values means the
# ---- Fourier series is not yet complete.

i = 0
while i < len(in_dihedral_coeffs):
    type_str = in_dihedral_coeffs[i][0]
    Kn = in_dihedral_coeffs[i][1]
    n = in_dihedral_coeffs[i][2]
    dn = in_dihedral_coeffs[i][3]

    #if (i>0):
    #    sys.stderr.write('prev_n='+str(in_dihedral_coeffs[i-1][-3])+'\n')
    #sys.stderr.write('n='+str(n)+'\n')

    if ((i>0) and (in_dihedral_coeffs[i-1][-3] < 0)):

        #sys.stdout.write('interaction_before_append: '+str(in_dihedral_coeffs[i-1])+'\n')
        assert(in_dihedral_coeffs[i-1][0] == in_dihedral_coeffs[i][0])
        in_dihedral_coeffs[i-1][-3] = -in_dihedral_coeffs[i-1][-3]
        comments = in_dihedral_coeffs[i-1][-1]
        in_dihedral_coeffs[i-1][-1] = Kn
        in_dihedral_coeffs[i-1].append(n)
        in_dihedral_coeffs[i-1].append(dn)
        in_dihedral_coeffs[i-1].append(comments)
        #sys.stdout.write('interaction_after_append: '+str(in_dihedral_coeffs[i-1])+'\n')
        del in_dihedral_coeffs[i]

    #elif len(in_dihedral_coeffs) < 3:
    #    del in_dihedral_coeffs[i]
    else:
        i += 1



for i in range(0, len(in_dihedral_coeffs)):
    type_str = in_dihedral_coeffs[i][0]
    params = in_dihedral_coeffs[i][1:]
    params = map(str, params)
    num_fourier_terms = (len(params)-1)/3
    dihedral_coeff_str = 'dihedral_coeff '+type_str+' '+\
        dihedral_style_name+' '+ \
        str(num_fourier_terms)+' '+ \
        ' '.join(params)
    in_dihedral_coeffs[i] = dihedral_coeff_str

# ---- finished processing dihedral fourier series ----


sys.stdout.write('  write_once(\"In Settings\") {\n    ')
sys.stdout.write('\n    '.join(in_dihedral_coeffs)+'\n')
sys.stdout.write('  } # (end of dihedral_coeffs)\n')





sys.stdout.write('\n')

sys.stdout.write('  write_once("Data Dihedrals By Type") {\n')

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
    dihedraltype = '@dihedral:'+atype1+'-'+atype2+'-'+atype3+'-'+atype4

    sys.stdout.write('    '+dihedraltype+' @atom:'+at1+' @atom:'+at2+' @atom:'+at3+' @atom:'+at4+'\n')

sys.stdout.write('  } # (end of Dihedrals By Type)\n')
sys.stdout.write('\n')


"""
         - 6 -      ***** INPUT FOR DIHEDRAL PARAMETERS *****

                    IPT , JPT , KPT , LPT , IDIVF , PK , PHASE , PN

                        FORMAT(A2,1X,A2,1X,A2,1X,A2,I4,3F15.2)

         IPT, ...   The atom symbols for the atoms forming a dihedral
                    angle.  If IPT .eq. 'X ' .and. LPT .eq. 'X ' then
                    any dihedrals in the system involving the atoms "JPT" and
                    and "KPT" are assigned the same parameters.  This is
                    called the general dihedral type and is of the form
                    "X "-"JPT"-"KPT"-"X ".

         IDIVF      The factor by which the torsional barrier is divided.
                    Consult Weiner, et al., JACS 106:765 (1984) p. 769 for
                    details. Basically, the actual torsional potential is

                           (PK/IDIVF) * (1 + cos(PN*phi - PHASE))

         PK         The barrier height divided by a factor of 2.

         PHASE      The phase shift angle in the torsional function.

                    The unit is degrees.

         PN         The periodicity of the torsional barrier.
                    NOTE: If PN .lt. 0.0 then the torsional potential
                          is assumed to have more than one term, and the
                          values of the rest of the terms are read from the
                          next cards until a positive PN is encountered.  The
                          negative value of pn is used only for identifying
                          the existence of the next term and only the
                          absolute value of PN is kept.

            The input is terminated by a blank card.
"""
