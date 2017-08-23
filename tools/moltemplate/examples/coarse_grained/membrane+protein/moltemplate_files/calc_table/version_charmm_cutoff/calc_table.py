#!/usr/bin/env python

# Calculate a table of pairwise energies and forces between "INT" atoms
# in the lipid membrane model described in
#   Brannigan et al, Phys Rev E, 72, 011915 (2005)
# The energy of this interaction U(r) = eps*(0.4*(sigma/r)^12 - 3.0*(sigma/r)^2)
# I realized later this is not what we want because although energy is conserved
# all enrgies are shifted with respect to energies used in the Brannigan paper
# (by 0.27 kCal/mole) and the later Watson JCP 2011 paper (by 0.224 kCal/mole).
# (So don't use this.)

# Calculate and print a

def S(r, rc1, rc2, derivative=False):
    """
    Calculate the switching function S(r) which decays continuously
    between 1 and 0 in the range from rc1 to rc2 (rc2>rc1):
       S(r) = (rc2^2 - r^2)^2 * (rc2^2 + 2*r^2 - 3*rc1^2) / (rc2^2-rc1^2)^3
    I'm using the same smoothing/switching cutoff function used by the CHARMM
    force-fields.  (I'm even using the same code to implement it, taken
    from lammps charmm/coul/charmm pair style, rewritten in python.)

    """
    assert(rc2>rc1)
    rsq           = r*r
    rc1sq         = rc1*rc1
    rc2sq         = rc2*rc2
    denom_lj_inv  = (1.0 / ((rc2sq-rc1sq)*
                            (rc2sq-rc1sq)*
                            (rc2sq-rc1sq)))
    if rsq > rc2sq:
        return 0.0
    elif rsq < rc1sq:
        if derivative:
            return 0.0
        else:
            return 1.0
    else:
        rc2sq_minus_rsq    = (rc2sq - rsq)
        rc2sq_minus_rsq_sq = rc2sq_minus_rsq * rc2sq_minus_rsq
        if derivative:
            return (12.0 * rsq * rc2sq_minus_rsq * (rsq-rc1sq) * denom_lj_inv)
        else:
            return (rc2sq_minus_rsq_sq *
                    (rc2sq + 2.0*rsq - 3.0*rc1sq) * denom_lj_inv)


def U(r, eps, sigma):
    return eps*   (0.4*pow((sigma/r),12) -   3.0*sigma*sigma/(r*r))

def F(r, eps, sigma):
    return eps*(12*0.4*pow((sigma/r),13)/sigma - 2*3.0*sigma*sigma/(r*r*r))

epsilon = 2.75/4.184 # kCal/mole
sigma   = 7.5
Rmin    = 0.02
Rmax    = 22.6
Rc1     = 22.0
Rc2     = 22.5
N       = 1130

for i in range(0,N):
    r = Rmin + i*(Rmax-Rmin)/(N-1)
    U_r = U(r, epsilon, sigma)
    F_r = F(r, epsilon, sigma)
    # Multiply U(r) & F(r) by the smoothing/switch function
    U_r = U_r * S(r, Rc1, Rc2)
    F_r = U_r * S(r, Rc1, Rc2, True)  +  F_r * S(r, Rc1, Rc2, False)
    print(str(i+1)+' '+str(r)+' '+str(U_r)+' '+str(F_r))

