#!/usr/bin/env python3
"""Generate a paired FCC system for fix symmetry under P-1 (inversion).

Writes two files:
  data.fcc_pinv     -- LAMMPS data file with explicit atom positions
  fcc_pinv.json     -- fix symmetry data file declaring the orbits

The FCC lattice is generated with the conventional 4-atom basis at
  (0,0,0), (0.5,0.5,0), (0.5,0,0.5), (0,0.5,0.5)
in units of the lattice parameter a, with the box covering nx*ny*nz
conventional cells starting at the origin. With this layout the lattice
is centrosymmetric through the cell origin (lamda = 0). Atoms whose
inversion image is themselves (mod 1) -- the eight half-integer-lamda
corner-sublattice atoms in a centered grid -- are written as Wyckoff
orbits; the remaining atoms come in inversion pairs.

The output of this script is what the in.lj_fcc_pinv example consumes.
Re-run it with different (nx, ny, nz) or lattice parameter to scale the
example to a different system size.
"""

import json
import sys

# --- Parameters (matches in.lj_fcc_pinv) ---------------------------------
nx, ny, nz = 4, 4, 4         # number of conventional FCC cells
lattice_a   = 1.5874         # = (4 / 0.8442)^(1/3) -- LJ density 0.8442 / atom
output_data = "data.fcc_pinv"
output_json = "fcc_pinv.json"

# --- Build atom list (basis-ordered to match LAMMPS create_atoms) --------
basis = [(0.0, 0.0, 0.0),
         (0.5, 0.5, 0.0),
         (0.5, 0.0, 0.5),
         (0.0, 0.5, 0.5)]
Lx, Ly, Lz = nx * lattice_a, ny * lattice_a, nz * lattice_a
atoms = []
for iz in range(nz):
    for iy in range(ny):
        for ix in range(nx):
            for bx, by, bz in basis:
                x = (ix + bx) * lattice_a
                y = (iy + by) * lattice_a
                z = (iz + bz) * lattice_a
                atoms.append((x, y, z))

natoms = len(atoms)
print(f"Generated {natoms} FCC atoms in [0,{Lx:.4f}]^3.")

# --- Pair atoms by inversion through cell origin (lamda = 0) -------------
# Lookup table from rounded fractional coords to tag.
def lamda(xyz):
    return (xyz[0] / Lx, xyz[1] / Ly, xyz[2] / Lz)

def key(s):
    return tuple(round(c % 1.0, 8) for c in s)

by_key = {}
for i, xyz in enumerate(atoms):
    by_key[key(lamda(xyz))] = i + 1   # LAMMPS tags are 1-based

orbits = []
paired = set()
n_wyckoff = 0
for tag in range(1, natoms + 1):
    if tag in paired:
        continue
    s = lamda(atoms[tag - 1])
    s_inv = tuple((-c) % 1.0 for c in s)
    if key(s) == key(s_inv):
        # Self-image: atom sits at an inversion center -> Wyckoff orbit.
        orbits.append({"tags": [tag, tag], "site_symmetry": [2]})
        paired.add(tag)
        n_wyckoff += 1
    else:
        partner = by_key[key(s_inv)]
        orbits.append({"tags": [tag, partner]})
        paired.update([tag, partner])

print(f"Orbits: {len(orbits)} total ({n_wyckoff} Wyckoff, "
      f"{len(orbits) - n_wyckoff} pairs).")

# --- Write LAMMPS data file ----------------------------------------------
with open(output_data, "w") as f:
    f.write("# FCC lattice, P-1 inversion through cell origin\n\n")
    f.write(f"{natoms} atoms\n1 atom types\n\n")
    f.write(f"0.0 {Lx:.10f} xlo xhi\n")
    f.write(f"0.0 {Ly:.10f} ylo yhi\n")
    f.write(f"0.0 {Lz:.10f} zlo zhi\n\n")
    f.write("Masses\n\n1 1.0\n\nAtoms # atomic\n\n")
    for i, (x, y, z) in enumerate(atoms):
        f.write(f"{i+1} 1 {x:.10f} {y:.10f} {z:.10f}\n")

# --- Write symmetry-data JSON --------------------------------------------
sym_data = {
    "group_name": "P-1",
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
