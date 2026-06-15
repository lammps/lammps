#!/usr/bin/env python3
"""Generate a LAMMPS data file of rigid 8-atom cubes on a 3D grid.

Each rigid body is a unit cube of 8 atoms (edge = `edge`). Bodies are placed on
a simple cubic grid with spacing `spacing` (> edge + cutoff) so that no two
bodies overlap or interact intramolecularly, which keeps the test well behaved
yet exercises body migration/exchange across processors.

Used to build a moderately sized test for fix rigid/small (and rigid/small/kk).
With n=5 this is 125 cubes / 1000 atoms (~70 KB), a reduced version of the
1000-cube system the PR author used to find migration/exchange bugs.
"""
import sys

n = int(sys.argv[1]) if len(sys.argv) > 1 else 5   # cubes per dimension
edge = 1.0          # cube edge length
spacing = 4.0       # grid spacing between cube centers (> edge + 2.5 cutoff/...)
mass = 1.0
out = sys.argv[2] if len(sys.argv) > 2 else "data.rigid.cubes"

# 8 corner offsets of a cube centered at the origin
h = edge / 2.0
corners = [(sx*h, sy*h, sz*h) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]

box = n * spacing
atoms = []
mol = 0
for i in range(n):
    for j in range(n):
        for k in range(n):
            mol += 1
            cx = (i + 0.5) * spacing
            cy = (j + 0.5) * spacing
            cz = (k + 0.5) * spacing
            for (ox, oy, oz) in corners:
                atoms.append((mol, 1, cx+ox, cy+oy, cz+oz))

with open(out, "w") as f:
    f.write("LAMMPS data file: %d rigid 8-atom cubes\n\n" % (n*n*n))
    f.write("%d atoms\n" % len(atoms))
    f.write("1 atom types\n\n")
    f.write("0.0 %g xlo xhi\n" % box)
    f.write("0.0 %g ylo yhi\n" % box)
    f.write("0.0 %g zlo zhi\n\n" % box)
    f.write("Masses\n\n1 %g\n\n" % mass)
    f.write("Atoms # molecular\n\n")
    for aid, (m, t, x, y, z) in enumerate(atoms, start=1):
        f.write("%d %d %d %.6f %.6f %.6f\n" % (aid, m, t, x, y, z))

print("wrote %s: %d atoms, %d bodies, box %.3f" % (out, len(atoms), mol, box))
