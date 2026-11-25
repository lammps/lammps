#!/usr/bin/env python3
"""Compare the thermo rows of two LAMMPS screen dumps.

A column is considered equal when |a-b| <= ATOL + RTOL*max(|a|,|b|).  The
absolute term matters: several decks print quantities that sit at round-off
around zero (a constraint force that should cancel, the fmax of a converged
minimization), where a pure relative test reports a huge difference for an
absolute difference of 1e-16.  The reported number is the ratio to that
tolerance, so 1.0 is exactly at the limit and anything above it fails.
"""
import sys, re, math

ATOL = 1.0e-8
RTOL = 1.0e-6

num = re.compile(r'^\s*\d+((\s+-?(\d+\.?\d*|nan|inf)([eE][-+]?\d+)?){2,})\s*$')

def rows(f):
    out = []
    try: lines = open(f, errors="replace")
    except OSError: return out
    for line in lines:
        if "%" in line or "CPU" in line: continue
        if num.match(line): out.append([float(x) for x in line.split()])
    return out

a, b, label = rows(sys.argv[1]), rows(sys.argv[2]), sys.argv[3]
if not a or not b:
    print("%-46s %-12s %s" % (label, "NO-OUTPUT", "cpu=%d kk=%d rows" % (len(a), len(b)))); sys.exit()
if len(a) != len(b):
    print("%-46s %-12s %s" % (label, "ROWMISMATCH", "cpu=%d kk=%d" % (len(a), len(b)))); sys.exit()
m = 0.0; nan = 0
for ra, rb in zip(a, b):
    for u, v in zip(ra, rb):
        if math.isnan(u) or math.isnan(v):
            if not (math.isnan(u) and math.isnan(v)): nan += 1
            continue
        m = max(m, abs(u - v) / (ATOL + RTOL * max(abs(u), abs(v))))
verdict = "OK" if (m <= 1.0 and nan == 0) else "DIVERGED"
print("%-46s %-12s %.3e%s" % (label, verdict, m, "  nan-mismatch=%d" % nan if nan else ""))
