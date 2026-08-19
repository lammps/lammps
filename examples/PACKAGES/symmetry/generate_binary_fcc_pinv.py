#!/usr/bin/env python3
"""Generate a binary FCC system for fix symmetry under P-1.

Writes:
  data.binary_fcc_pinv  -- LAMMPS data file with two atom types
  binary_fcc_pinv.json  -- fix symmetry data file declaring the orbits

This mirrors generate_fcc_pinv.py but assigns two atom types in an L1_2-style
ordered pattern: the corner-sublattice atoms are type 1 (one per conventional
cell) and the three face-center sublattices are type 2 (three per cell), so
the type ratio is 1:3. Inversion through the cell origin pairs each atom
with one of the same type (FCC sublattices are preserved by inversion), so
the orbits split cleanly into "type-1 orbits" and "type-2 orbits".
"""

import json

# --- Parameters ----------------------------------------------------------
N           = 4              # cells per dimension (nx = ny = nz = N)
lattice_a   = 1.5874         # same as in.lj_fcc_pinv (LJ density ~0.85)
output_data = "data.binary_fcc_pinv"
output_json = "binary_fcc_pinv.json"

# --- Build atom list: type 1 = corner sublattice, type 2 = face-centers --
# Basis ordering is preserved across cells so atom tags are predictable.
basis = [
    ((0.00, 0.00, 0.00), 1),
    ((0.50, 0.50, 0.00), 2),
    ((0.50, 0.00, 0.50), 2),
    ((0.00, 0.50, 0.50), 2),
]
L = N * lattice_a
atoms = []      # list of ((x, y, z), type)
for iz in range(N):
    for iy in range(N):
        for ix in range(N):
            for (bx, by, bz), t in basis:
                atoms.append((((ix + bx) * lattice_a,
                               (iy + by) * lattice_a,
                               (iz + bz) * lattice_a), t))
natoms = len(atoms)
ntype1 = sum(1 for _, t in atoms if t == 1)
ntype2 = sum(1 for _, t in atoms if t == 2)
print(f"Generated {natoms} atoms in [0,{L:.4f}]^3: "
      f"{ntype1} type-1, {ntype2} type-2.")

# --- Pair atoms by inversion through cell origin (lamda = 0) -------------
def lamda(xyz):
    return (xyz[0] / L, xyz[1] / L, xyz[2] / L)

def key(s):
    return tuple(round(c % 1.0, 8) for c in s)

by_key = {key(lamda(xyz)): (i + 1, t) for i, (xyz, t) in enumerate(atoms)}

orbits = []
paired = set()
n_wyckoff = 0
for tag in range(1, natoms + 1):
    if tag in paired:
        continue
    xyz, t = atoms[tag - 1]
    s = lamda(xyz)
    s_inv = tuple((-c) % 1.0 for c in s)
    partner, partner_type = by_key[key(s_inv)]
    if partner_type != t:
        raise SystemExit(
            f"Inversion would map atom {tag} (type {t}) to atom {partner} "
            f"(type {partner_type}). The structure is not inversion-symmetric "
            "in this type assignment."
        )
    if key(s) == key(s_inv):
        orbits.append({"tags": [tag, tag], "site_symmetry": [2]})
        paired.add(tag)
        n_wyckoff += 1
    else:
        orbits.append({"tags": [tag, partner]})
        paired.update([tag, partner])

print(f"Orbits: {len(orbits)} total "
      f"({n_wyckoff} Wyckoff, {len(orbits) - n_wyckoff} pairs).")

# --- Write LAMMPS data file ----------------------------------------------
with open(output_data, "w") as f:
    f.write("# Binary FCC (1:3 ratio), P-1 inversion through cell origin\n\n")
    f.write(f"{natoms} atoms\n2 atom types\n\n")
    f.write(f"0.0 {L:.10f} xlo xhi\n")
    f.write(f"0.0 {L:.10f} ylo yhi\n")
    f.write(f"0.0 {L:.10f} zlo zhi\n\n")
    f.write("Masses\n\n1 1.0\n2 1.0\n\nAtoms # atomic\n\n")
    for i, ((x, y, z), t) in enumerate(atoms):
        f.write(f"{i+1} {t} {x:.10f} {y:.10f} {z:.10f}\n")

# --- Write symmetry-data JSON --------------------------------------------
sym_data = {
    "group_name": "P-1 (binary L1_2-like)",
    "lattice":    "triclinic",
    "ops": [
        {"R": [[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]], "t": [0, 0, 0]},
        {"R": [[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]], "t": [0, 0, 0]},
    ],
    "orbits": orbits,
}
with open(output_json, "w") as f:
    json.dump(sym_data, f, indent=2)
    f.write("\n")

print(f"Wrote {output_data} and {output_json}.")
