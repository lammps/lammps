#!/usr/bin/env python3
"""Generate a SiC zinc-blende system constrained by a 2-fold rotation about z.

Writes:
  data.sic_c2   -- LAMMPS data file (atom type 1 = C, type 2 = Si)
  sic_c2.json   -- fix symmetry data file declaring the orbits

Zinc-blende SiC is *not* centrosymmetric -- inversion would swap Si and C
because the two FCC sublattices sit at distinct positions
((0,0,0) for Si vs (1/4,1/4,1/4) for C). P-1 inversion therefore cannot
be used here. Instead we pick a single type-preserving element of the
F-43m point group: the 2-fold rotation about the cubic z-axis,
(x, y, z) -> (-x, -y, z), which keeps each FCC sublattice on itself.

  Op 1: identity
  Op 2: R = diag(-1, -1, 1), t = 0

Atoms whose (lamda_x, lamda_y) both land in {0, 1/2} are invariant under
op 2 -- those become Wyckoff orbits whose asym atom is constrained to
the rotation axis (x and y pinned, z free). In a 2x2x2 conventional cell
that's the 8 Si atoms on the corner sublattice; everything else comes in
2-fold pairs. The C sublattice (offset by (1/4, 1/4, 1/4) from Si) has
no Wyckoff atoms under this rotation -- all C atoms come in pairs.

The atom-type ordering (1 = C, 2 = Si) matches the convention of the
SiC.tersoff parameter file used by in.sic_c2.
"""

import json

# --- Parameters ----------------------------------------------------------
N           = 2              # cells per dimension (nx = ny = nz = N)
lattice_a   = 4.360          # SiC conventional lattice parameter (Angstrom)
output_data = "data.sic_c2"
output_json = "sic_c2.json"

# --- Build atom list with two types (1 = C, 2 = Si) ----------------------
# Si on the FCC sublattice (type 2); C on the FCC + (1/4, 1/4, 1/4) sublattice
# (type 1). Atoms ordered by (cell, basis) so tags match standard create_atoms.
basis = [
    ((0.00, 0.00, 0.00), 2),    # Si
    ((0.00, 0.50, 0.50), 2),    # Si
    ((0.50, 0.00, 0.50), 2),    # Si
    ((0.50, 0.50, 0.00), 2),    # Si
    ((0.25, 0.25, 0.25), 1),    # C
    ((0.25, 0.75, 0.75), 1),    # C
    ((0.75, 0.25, 0.75), 1),    # C
    ((0.75, 0.75, 0.25), 1),    # C
]
L = N * lattice_a
atoms = []
for iz in range(N):
    for iy in range(N):
        for ix in range(N):
            for (bx, by, bz), t in basis:
                atoms.append((((ix + bx) * lattice_a,
                               (iy + by) * lattice_a,
                               (iz + bz) * lattice_a), t))
natoms = len(atoms)
n_si = sum(1 for _, t in atoms if t == 2)
n_c  = sum(1 for _, t in atoms if t == 1)
print(f"Generated {natoms} SiC atoms in [0,{L:.4f}]^3: "
      f"{n_si} Si (type 2), {n_c} C (type 1).")

# --- Pair atoms by 2-fold rotation around z: s' = (-s_x, -s_y, s_z) ------
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
    s_rot = ((-s[0]) % 1.0, (-s[1]) % 1.0, s[2] % 1.0)
    partner, partner_type = by_key[key(s_rot)]
    if partner_type != t:
        raise SystemExit(
            f"Op would map atom {tag} (type {t}) to atom {partner} "
            f"(type {partner_type}) -- rotation is not type-preserving here."
        )
    if key(s) == key(s_rot):
        orbits.append({"tags": [tag, tag], "site_symmetry": [2]})
        paired.add(tag)
        n_wyckoff += 1
    else:
        orbits.append({"tags": [tag, partner]})
        paired.update([tag, partner])

print(f"Orbits: {len(orbits)} total "
      f"({n_wyckoff} Wyckoff on the rotation axis, "
      f"{len(orbits) - n_wyckoff} 2-fold pairs).")

# --- Write LAMMPS data file ----------------------------------------------
with open(output_data, "w") as f:
    f.write("# SiC zinc-blende, C_2 rotation about z\n")
    f.write("# atom type 1 = C, atom type 2 = Si\n\n")
    f.write(f"{natoms} atoms\n2 atom types\n\n")
    f.write(f"0.0 {L:.10f} xlo xhi\n")
    f.write(f"0.0 {L:.10f} ylo yhi\n")
    f.write(f"0.0 {L:.10f} zlo zhi\n\n")
    f.write("Masses\n\n1 12.0107\n2 28.0855\n\nAtoms # atomic\n\n")
    for i, ((x, y, z), t) in enumerate(atoms):
        f.write(f"{i+1} {t} {x:.10f} {y:.10f} {z:.10f}\n")

# --- Write symmetry-data JSON --------------------------------------------
sym_data = {
    "group_name": "C_2(z) subgroup of F-43m (SiC zinc-blende)",
    "lattice":    "cubic",
    "ops": [
        {"R": [[ 1,  0, 0], [ 0,  1, 0], [0, 0, 1]], "t": [0, 0, 0]},
        {"R": [[-1,  0, 0], [ 0, -1, 0], [0, 0, 1]], "t": [0, 0, 0]},
    ],
    "orbits": orbits,
}
with open(output_json, "w") as f:
    json.dump(sym_data, f, indent=2)
    f.write("\n")

print(f"Wrote {output_data} and {output_json}.")
