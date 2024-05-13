#!/usr/bin/env python

from tabulate import PairTabulate

################################################################################
import math

# this table contains a correction to be added to the Al_zhou.eam.alloy potential via hybrid/overlay
# the purpose is to smoothly replace the original pairwise repulsion with the ZBL potential.
# the combined potential switches from ZBL to EAM w/o embedding and then subtracts the full EAM term.
# this way the table can be added to the eam/alloy pair style via hybrid/overlay
# due to using a Fermi-like switching function (tanh()) there are no discontinuities in energy or force

def eam_2body(r):
    biga  = 0.251519 # metal units
    bigb  = 0.313394 # metal units
    rzero = 2.886166
    alpha = 6.942419
    beta  = 3.702623
    lamda = 0.790264
    kappa = 0.395132
    return (biga * math.exp(-alpha * (r/rzero - 1.0))) / (1.0 + math.pow(r/rzero - kappa, 20)) \
         - (bigb * math.exp(-beta  * (r/rzero - 1.0))) / (1.0 + math.pow(r/rzero - lamda, 20))

def zbl_energy(r):
    qqr2e = 14.399645  # for metal units
    pzbl = 0.23
    a0 = 0.46850
    c1 = 0.02817
    c2 = 0.28022
    c3 = 0.50986
    c4 = 0.18175
    d1 = 0.20162
    d2 = 0.40290
    d3 = 0.94229
    d4 = 3.19980

    zi = 13.0  # aluminium
    zj = 13.0  # aluminium

    prefactor = zi * zj * qqr2e / r
    rbya = r * (math.pow(zi, pzbl) + math.pow(zj, pzbl)) / a0
    f = prefactor * (c4*math.exp(-d4*rbya) + c3*math.exp(-d3*rbya) \
                     + c2*math.exp(-d2*rbya) + c1*math.exp(-d1*rbya))
    return f

def combined(r):
    rmid = 2.886166
    switch_on = 0.5*(math.tanh(math.exp(1.0)*(r - rmid))+1.0)
    return (1.0-switch_on)*zbl_energy(r) + (switch_on-1.0)*eam_2body(r)

################################################################################

if __name__ == "__main__":
    ptable = PairTabulate(combined, units='metal', \
                          comment='Correct Al_zhou EAM with ZBL using hybrid/overlay')
    ptable.run('AL_EAM_ZBL_CORR')

# NOTE use/create table with an outer cutoff of 10.1025 to match the EAM potential cutoff.
