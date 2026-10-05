#!/usr/bin/env python
# Convert the output of a LAMMPS 'fix ave/time ... mode vector' applied to
# 'compute mbar' into a u_kln array (nstates, nstates, nsamples) for pymbar.
#
# compute mbar writes, at every sampled step, the reduced potentials U_l/kT of
# the current configuration evaluated at every state l. With a ramped fix
# adapt/fep the configuration sampled at step t belongs to the held state
# k = (t - 1) // window. Grouping the samples by k and stacking the per-step
# vectors yields u_kln[k, l, n] = (sample n drawn from state k) evaluated at l,
# which is what mbar.py expects.

import sys
from argparse import ArgumentParser
import numpy as np


def sigdigits(text):
    """number of significant digits in the text of a floating-point number"""
    mantissa = text.lower().split('e')[0].lstrip('+-').replace('.', '')
    return len(mantissa.lstrip('0'))


parser = ArgumentParser(description='Reshape LAMMPS compute mbar output into u_kln for pymbar.')
parser.add_argument('infile', help='LAMMPS fix ave/time vector file (e.g. mbar.lmp)')
parser.add_argument('window', type=int, help='length in steps of each held-lambda stage (e.g. 50000)')
parser.add_argument('outfile', help='output .npy file with the u_kln array')
args = parser.parse_args()

# parse blocks: each output is a header line "timestep nrows" followed by
# nrows lines "row value"; comment lines (#) appear only in the file header
samples = {}        # state index k -> list of per-step reduced-potential vectors
nstates = None
ndigits = 0         # largest number of significant digits found in the values
with open(args.infile) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        step, nrows = line.split()
        step, nrows = int(step), int(nrows)
        if nstates is None:
            nstates = nrows
        vec = np.empty(nrows)
        for _ in range(nrows):
            row = next(f).split()
            if len(row) != 2:
                raise SystemExit(f'malformed line "{" ".join(row)}" at step {step}: the format '
                                 'of fix ave/time must start with a space, e.g. " %.15g"')
            idx, val = row
            vec[int(idx) - 1] = float(val)
            ndigits = max(ndigits, sigdigits(val))
        # state held during step t (fix adapt/fep updates after each window);
        # step 0 is the initial equilibrated config, which belongs to state 0
        k = max(0, (step - 1) // args.window)
        samples.setdefault(k, []).append(vec)

ks = sorted(samples)
if ks != list(range(nstates)):
    raise SystemExit(f'expected states 0..{nstates-1}, found {ks} '
                     f'(check that window={args.window} matches the stage length)')

nsamp = min(len(samples[k]) for k in ks)       # truncate to equal counts
u_kln = np.zeros((nstates, nstates, nsamp))
for k in ks:
    u_kln[k, :, :] = np.array(samples[k][:nsamp]).T   # (nstates, nsamp)

# the reduced potentials are large numbers and MBAR depends on their small
# differences, so the 6 digits of the default format of fix ave/time are not enough
if ndigits <= 6:
    print(f'WARNING: the reduced potentials in {args.infile} have at most {ndigits} '
          'significant digits, which degrades the accuracy of MBAR. Write them with '
          'more digits, e.g. with format " %.15g" in fix ave/time.', file=sys.stderr)

np.save(args.outfile, u_kln)
print(f'states = {nstates}, samples/state = {nsamp} (min), saved {args.outfile}')
