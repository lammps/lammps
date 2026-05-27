#!/usr/bin/env python3
"""Generate a paired diamond-lattice system for fix symmetry under P-1.

Writes:
  data.diamond_pinv  -- LAMMPS data file with explicit atom positions
  diamond_pinv.json  -- fix symmetry data file declaring the orbits

The diamond structure has 8 atoms per conventional cell, two interpenetrating
FCC sublattices offset by (1/4, 1/4, 1/4). It is centrosymmetric, but the
inversion center is *not* on an atom: it sits at (a/8, a/8, a/8) and
equivalent positions, where a is the conventional lattice parameter.

For nx*ny*nz conventional cells with box edge N*a (here N = nx = ny = nz),
inversion through Cartesian (a/8, a/8, a/8) corresponds to the affine
operator in fractional coordinates s' = R * s + t with

    R = -I,    t = (1/(4N), 1/(4N), 1/(4N))

Since no atom sits on an inversion center, every orbit is a clean pair --
no Wyckoff atoms.
"""

import json

# --- Parameters ----------------------------------------------------------
N           = 2              # cells per dimension (nx = ny = nz = N)
# Conventional diamond lattice parameter (LJ units). Diamond NN distance
# = a * sqrt(3) / 4; set a so NN sits near the LJ minimum 2^(1/6) sigma
# (with sigma = 1). a = 2.6 -> NN ~ 1.126.
lattice_a   = 2.6
output_data = "data.diamond_pinv"
output_json = "diamond_pinv.json"

# --- Build atom list (basis-ordered to match LAMMPS create_atoms) --------
# Diamond has 8 atoms per conventional cell: FCC + FCC shifted by (1/4,1/4,1/4).
basis = [(0.00, 0.00, 0.00),
         (0.00, 0.50, 0.50),
         (0.50, 0.00, 0.50),
         (0.50, 0.50, 0.00),
         (0.25, 0.25, 0.25),
         (0.25, 0.75, 0.75),
         (0.75, 0.25, 0.75),
         (0.75, 0.75, 0.25)]
L = N * lattice_a
atoms = []
for iz in range(N):
    for iy in range(N):
        for ix in range(N):
            for bx, by, bz in basis:
                atoms.append(((ix + bx) * lattice_a,
                              (iy + by) * lattice_a,
                              (iz + bz) * lattice_a))
natoms = len(atoms)
print(f"Generated {natoms} diamond atoms in [0,{L:.4f}]^3.")

# --- Pair atoms by inversion through (a/8, a/8, a/8) Cartesian -----------
t_frac = (1.0 / (4 * N), 1.0 / (4 * N), 1.0 / (4 * N))

def lamda(xyz):
    return (xyz[0] / L, xyz[1] / L, xyz[2] / L)

def key(s):
    return tuple(round(c % 1.0, 8) for c in s)

by_key = {key(lamda(xyz)): i + 1 for i, xyz in enumerate(atoms)}

orbits = []
paired = set()
for tag in range(1, natoms + 1):
    if tag in paired:
        continue
    s = lamda(atoms[tag - 1])
    s_inv = tuple((-s[k] + t_frac[k]) % 1.0 for k in range(3))
    if key(s) == key(s_inv):
        # would be Wyckoff, but the inversion center is off-atom in diamond
        # so this branch shouldn't fire -- guard it anyway.
        orbits.append({"tags": [tag, tag], "site_symmetry": [2]})
        paired.add(tag)
    else:
        partner = by_key.get(key(s_inv))
        if partner is None:
            raise SystemExit(
                f"Atom {tag} at lamda {s} has inversion image at {s_inv} "
                "but no matching atom found -- check N or the t vector"
            )
        orbits.append({"tags": [tag, partner]})
        paired.update([tag, partner])

print(f"Orbits: {len(orbits)} pairs (no Wyckoff -- diamond inversion centers "
      f"lie off-atom).")

# --- Write LAMMPS data file ----------------------------------------------
with open(output_data, "w") as f:
    f.write("# Diamond lattice, P-1 inversion through Cartesian (a/8, a/8, a/8)\n\n")
    f.write(f"{natoms} atoms\n1 atom types\n\n")
    f.write(f"0.0 {L:.10f} xlo xhi\n")
    f.write(f"0.0 {L:.10f} ylo yhi\n")
    f.write(f"0.0 {L:.10f} zlo zhi\n\n")
    f.write("Masses\n\n1 1.0\n\nAtoms # atomic\n\n")
    for i, (x, y, z) in enumerate(atoms):
        f.write(f"{i+1} 1 {x:.10f} {y:.10f} {z:.10f}\n")

# --- Write symmetry-data JSON --------------------------------------------
def fmt(x):
    """Emit small fractions as rational strings, else as numbers."""
    from fractions import Fraction
    fr = Fraction(x).limit_denominator(64)
    if abs(float(fr) - x) < 1e-12:
        return str(fr) if fr.denominator != 1 else fr.numerator
    return x

sym_data = {
    "group_name": "P-1 (diamond)",
    "lattice":    "triclinic",
    "ops": [
        {"R": [[ 1, 0, 0], [ 0, 1, 0], [ 0, 0, 1]], "t": [0, 0, 0]},
        {"R": [[-1, 0, 0], [ 0,-1, 0], [ 0, 0,-1]],
         "t": [fmt(t_frac[0]), fmt(t_frac[1]), fmt(t_frac[2])]},
    ],
    "orbits": orbits,
}
with open(output_json, "w") as f:
    json.dump(sym_data, f, indent=2)
    f.write("\n")

print(f"Wrote {output_data} and {output_json}.")
