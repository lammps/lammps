fix symmetry examples
=====================

These examples show how `fix symmetry` keeps a crystal, or any periodic
system, exactly symmetric under a chosen space group during MD and
minimization.  The method follows Cox and White, J. Chem. Theory
Comput. 18, 4077 (2022).  See the fix symmetry documentation in the
manual for details.

Each input builds its system with `create_atoms` and uses a small JSON
file that lists the symmetry operators.  The fix finds the sets of
symmetry related atoms (orbits) from the positions, so the file does not
contain atom IDs and works for any number of atoms.

Inputs
------

| Input                   | Atoms | Potential | Symmetry                                  | Shows |
|-------------------------|-------|-----------|-------------------------------------------|-------|
| `in.lj_pinv`            |     8 | LJ        | P-1 (inversion), `p-1.json`               | Why the fix is needed: without it, a tiny perturbation destroys the symmetry. |
| `in.lj_pm_wyckoff`      |     3 | LJ        | Pm (mirror plane), `pm.json`              | An atom on the mirror plane, a special position, moves only within the plane. |
| `in.lj_pmm2_wyckoff`    |     6 | LJ        | Pmm2, `pmm2.json`                         | Two atoms on a mirror plane form an orbit of 2 atoms for 4 operators. |
| `in.lj_fcc_pinv`        |   256 | LJ        | P-1, `p-1.json`                           | Thermostat with the correct number of degrees of freedom; 8 atoms on inversion centers. |
| `in.lj_binary_fcc_pinv` |   256 | LJ        | P-1, `p-1.json`                           | Two atom types (L1_2 structure); each orbit contains only one type. |
| `in.lj_diamond_pinv`    |    64 | LJ        | P-1 with a shifted center, `p-1_shifted.json` | An operator with a translation: inversion through a point between atoms. |
| `in.sic_c2`             |    64 | Tersoff   | P2 (two-fold axis), `p2.json`             | Zinc blende SiC, where inversion would exchange Si and C. |
| `in.lj_fcc_relax`       |   256 | LJ        | Fm-3m with 192 operators, `fm-3m.json`    | Box relaxation under the full cubic symmetry; file written by the find-symmetry tool with a supercell entry. |

Most inputs print a symmetry indicator in the last thermo column: the sum
of sin(2 pi s) over the fractional coordinates s of all atoms.  For each
atom and its image the terms cancel, so the sum is zero (to round-off)
for a symmetric configuration and grows when the symmetry is broken.

Symmetry files
--------------

A symmetry file lists the operators in the coordinate triplet notation
of the International Tables for Crystallography and of CIF files:

```json
{
  "group_name": "Pmm2",
  "ops": [ "x,y,z", "-x,-y,z", "-x,y,z", "x,-y,z" ]
}
```

The operators act on fractional coordinates of the simulation box and
must form a group, including the identity.  Entries starting with an
underscore are comments.  If the operators are those of a unit cell and
the box contains several unit cells, add `"supercell": [nx, ny, nz]`;
the fix then adds the translations between unit cells.

The `find-symmetry` tool in `tools/symmetry` writes such files: either
for the space group of the atoms in a LAMMPS data file, or for a space
group given by its number.  `fm-3m.json` was written with:

```
find-symmetry -g 225 -s 4 4 4 -o fm-3m.json
```

Setting up a run
----------------

* The initial positions must be symmetric within the tolerance of the
  fix (keyword `tol`, default 1.0e-6 in fractional coordinates).  The fix
  removes any remaining deviation at the first step.
* Random velocities from the `velocity` command are not symmetric.  A
  `run 0` with the fix makes them symmetric, which removes about half of
  the kinetic energy; `velocity all scale` restores the temperature.
  The temperature is computed from the degrees of freedom left by the
  symmetry, so thermostats work as usual.
* Define fix symmetry after all fixes that add or change forces.
* The box must fit the operators: for example a four-fold axis along z
  requires equal box lengths along x and y.  This is checked at the
  start and whenever the box changes.

Common pitfalls
---------------

* **Operators that exchange atom types.**  Each atom must have an image
  of the same type under every operator.  Zinc blende (SiC, GaAs, ...)
  has no inversion symmetry because inversion exchanges the two
  sublattices; `in.sic_c2` uses a two-fold axis instead.

* **Too much symmetry.**  Atoms on special positions with a site
  symmetry that fixes all coordinates cannot move.  With the full space
  group of a crystal and the translations between all unit cells, as in
  `in.lj_fcc_relax`, every atom is such a case and only the box can
  change.  Use a subgroup to leave more coordinates free.

* **Screw axes and glide planes in a box of several unit cells.**  Their
  operators combine to translations by a unit cell, so these translations
  must be part of the group: use the `supercell` entry.

* **Box shape.**  Hexagonal and trigonal operators such as `-y,x-y,z`
  need a triclinic box with equal lengths a and b at 120 degrees.

References
----------

* S. Cox and A. D. White, *Symmetric Molecular Dynamics*, J. Chem.
  Theory Comput. **18**, 4077-4081 (2022), DOI 10.1021/acs.jctc.2c00401
* Reference implementation symd: https://github.com/whitead/symd
