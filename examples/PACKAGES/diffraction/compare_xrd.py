#!/usr/bin/env python3
"""Compare two global arrays written by fix ave/time for compute xrd and
compute xrd/fft.

Usage:  python3 compare_xrd.py xrd_dft.out xrd_fft.out [tolerance]

The two computes must produce the same rows in the same order, so the 2theta
column has to match exactly.  The intensity column is compared both in absolute
terms scaled by the strongest reflection (the dynamic range that matters for a
diffraction pattern) and row by row in relative terms, with the weak rows
reported separately because they are where the spreading error shows up first.

Exits nonzero if the rows disagree or a tolerance is given and exceeded.
"""

import sys


def load(filename):
    """Return the (2theta, intensity) rows of a fix ave/time vector file."""
    rows = []
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split()
            # a blank line, or the "timestep number-of-rows" record
            if len(parts) < 3:
                continue
            rows.append((float(parts[1]), float(parts[2])))
    return rows


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1

    ref = load(sys.argv[1])
    new = load(sys.argv[2])
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else None

    if len(ref) != len(new):
        print("FAIL: %s has %d rows, %s has %d"
              % (sys.argv[1], len(ref), sys.argv[2], len(new)))
        return 1

    imax = max(row[1] for row in ref)
    if imax <= 0.0:
        print("FAIL: reference intensities are all zero")
        return 1

    maxabs = 0.0
    maxrel = 0.0
    maxrel_weak = 0.0
    nweak = 0
    worst = None

    for n, ((t_ref, i_ref), (t_new, i_new)) in enumerate(zip(ref, new)):
        if abs(t_ref - t_new) > 1.0e-12 * max(1.0, abs(t_ref)):
            print("FAIL: row %d has 2theta %.17g vs %.17g" % (n, t_ref, t_new))
            return 1

        diff = abs(i_ref - i_new)
        maxabs = max(maxabs, diff)

        # a relative error is only meaningful where the reference is above the
        # roundoff floor of the direct sum itself
        if i_ref > 1.0e-10 * imax:
            rel = diff / i_ref
            if rel > maxrel:
                maxrel = rel
                worst = (n, t_ref, i_ref, i_new)
            if i_ref < 1.0e-3 * imax:
                nweak += 1
                maxrel_weak = max(maxrel_weak, rel)

    print("rows compared            : %d" % len(ref))
    print("strongest reflection     : %.6g" % imax)
    print("max |dI| / I_max         : %.3g" % (maxabs / imax))
    print("max relative error       : %.3g" % maxrel)
    if worst:
        print("  worst row %d: 2theta %.6g, %.8g vs %.8g"
              % (worst[0], worst[1], worst[2], worst[3]))
    print("max relative error on the %d weak rows (I < 1e-3 I_max): %.3g"
          % (nweak, maxrel_weak))

    if tol is not None and maxabs / imax > tol:
        print("FAIL: max |dI| / I_max exceeds tolerance %.3g" % tol)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
