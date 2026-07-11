#!/usr/env/python3

import os.path as op
import sys


def rel_error(out, ref):
    return 100 * abs((ref - out) / out)


assert len(sys.argv) == 5

with open(sys.argv[1]) as f:
    ref_line = f.readlines()[-1].split(", ")
e_ref = float(ref_line[1])
q_ref = float(ref_line[3])
inv11_ref = float(ref_line[4])
inv12_ref = float(ref_line[5])
b1_ref = float(ref_line[6])
felec1_ref = float(ref_line[8])
felyt1_ref = float(ref_line[10])
press_ref = float(ref_line[12])

# out.csv
with open(sys.argv[2]) as f:
    out_line = f.readlines()[-1].split(", ")
e_out = float(out_line[0])
q_out = float(out_line[1])
press_out = float(out_line[2])

out_lines = [
    ("energy", e_ref, e_out),
    ("charge", q_ref, q_out),
    ("pressure", press_ref, press_out),
]

# vec.csv
vec_file = "vec.csv"
if op.isfile(vec_file):
    with open(sys.argv[4]) as f:
        vec_line = f.readlines()[1:]
    b1_out = float(vec_line[0])
    b2_out = float(vec_line[1])
    out_lines.append(("b1", b1_ref, b1_out))

# inv.csv
inv_file = "inv.csv"
if op.isfile(inv_file):
    with open(inv_file) as f:
        inv_line = f.readlines()[1].split()
    inv11_out = float(inv_line[0])
    inv12_out = float(inv_line[1])
    out_lines.append(("inv11", inv11_ref, inv11_out))

# forces.lammpstrj
force_file = "forces.lammpstrj"
with open(force_file) as f:
    lines = f.readlines()[9:]
    for name, i, f_ref in [("felec1", "1", felec1_ref), ("felyt1", "3", felyt1_ref)]:
        f_out = next(float(y[3]) for x in lines if (y := x.split()) and y[0] == i)
        out_lines.append((name, f_ref, f_out))

lines = []
for label, ref, out in out_lines:
    error = rel_error(out, ref)
    lines.append(f"{label}: {out:.5f}, {error:.5f}\n")

with open("madelung.txt", "a") as f:
    f.writelines(lines)
