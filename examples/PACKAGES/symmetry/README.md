fix symmetry examples
=====================

These examples demonstrate the use of `fix symmetry` on several
crystalline systems of increasing complexity. The fix constrains a
molecular-dynamics run to obey a chosen crystallographic space-group
symmetry at every timestep, following Cox & White, J. Chem. Theor.
Comput. 18, 4077 (2022). See the fix symmetry documentation in the
manual for further details.

Provided examples
-----------------

In rough order of complexity. Each input prints a residual check at the
end of the run; the constraint invariants must stay at the FP noise
floor (typically 1e-16 to 1e-12, depending on the integrator length).

| Input                   | Atoms | Pair style | Symmetry op(s)                 | Notes |
|-------------------------|-------|------------|--------------------------------|-------|
| `in.lj_pinv`            |     8 | lj/cut     | P-1 (inversion through origin) | Smallest pedagogical example. Four asymmetric atoms in one half of the box, four inversion images in the other half. |
| `in.lj_pm_wyckoff`      |     3 | lj/cut     | Pm (identity + y-mirror)       | Two atoms in a general-position mirror pair plus one atom on the mirror plane (Wyckoff site). |
| `in.lj_fcc_pinv`        |   256 | lj/cut     | P-1                            | 4x4x4 FCC crystal; **132 orbits** with 124 inversion pairs + 8 Wyckoff atoms at the corner-sublattice half-integer-lamda positions. |
| `in.lj_diamond_pinv`    |    64 | lj/cut     | P-1 with off-origin inversion  | 2x2x2 diamond crystal. Inversion center sits at Cartesian `(a/8, a/8, a/8)` -- *between* atoms, so every orbit is a clean pair, no Wyckoff. |
| `in.lj_binary_fcc_pinv` |   256 | lj/cut     | P-1                            | Same FCC layout as `in.lj_fcc_pinv` but with a 1:3 L1_2-style ordered binary alloy. Three `pair_coeff` lines: stronger 1-1, weaker 2-2, mixed 1-2. FCC sublattices are inversion-invariant, so orbits are type-pure. |
| `in.sic_c2`             | 64    | tersoff    | C_2(z) subgroup of F-43m       | 2x2x2 zinc-blende SiC. Inversion would swap Si and C, so P-1 cannot be used here; the 2-fold rotation about z preserves both sublattices. 8 Si atoms sit on the rotation axis (Wyckoff: constrained to move only along z); 28 orbits are 2-fold pairs across z. |

Generator scripts
-----------------

Each non-trivial example ships a Python script that emits its data file
and JSON in one shot. Re-run after editing the parameters at the top of
the script to scale to a different system size or lattice parameter.

| Script                        | Produces                                       | Settings at top              |
|-------------------------------|------------------------------------------------|------------------------------|
| `generate_fcc_pinv.py`        | `data.fcc_pinv`, `fcc_pinv.json`               | `N` (cells/dim), `lattice_a` |
| `generate_diamond_pinv.py`    | `data.diamond_pinv`, `diamond_pinv.json`       | `N`, `lattice_a`             |
| `generate_binary_fcc_pinv.py` | `data.binary_fcc_pinv`, `binary_fcc_pinv.json` | `N`, `lattice_a`             |
| `generate_sic_c2.py`          | `data.sic_c2`, `sic_c2.json`                   | `N`, `lattice_a`             |

All scripts follow the same three-phase structure:

1. **Build atom list** -- iterate cells × basis and record `((x, y, z), type)`.
2. **Pair atoms by the chosen symmetry op** -- look up each atom's image
   in a fractional-coordinate hash table; flag self-images as Wyckoff;
   error out if the op would swap atom types.
3. **Write data file + JSON** -- the JSON encodes the operator(s) and
   the orbit table (`tags` arrays parallel to `ops`, plus optional
   `site_symmetry` for Wyckoff orbits).

The standalone `symgroup.p1.json` and `symgroup.pm.json` are operator-only
stubs (no orbits) intended as starting points for hand-edited files.

Adapting these examples to your own crystal
-------------------------------------------

Three pieces typically need to change. Working from a generator script
closest to your system is usually easier than starting from scratch.

1. **Atom layout (basis + cells)** -- edit the `basis` list at the top
   of the generator. Each entry is `((x_frac, y_frac, z_frac), type)`
   in conventional-cell fractional coordinates. The `N` and `lattice_a`
   constants then scale the structure.

2. **Symmetry operator(s)** -- replace the `s_inv = ...` line that
   computes each atom's image under your chosen operation, and update
   the `"ops"` block in the JSON-emit section accordingly. Reminders:

   * `R` is a 3x3 integer matrix with `|det R| = 1`.
   * `t` is a length-3 fractional translation; numbers or rational
     strings (`"1/2"`) are accepted by the parser.
   * `ops[0]` (the first operator) must be the identity; the fix
     verifies group closure under composition modulo lattice translations.
   * The operator must be **type-preserving** for multi-type systems --
     the generators include an explicit check that errors out if an
     atom's image is of a different type.

3. **Orbit and Wyckoff structure** -- the generators auto-detect
   Wyckoff sites (atoms that map onto themselves under the operator)
   and write `"site_symmetry": [op_index]` on those orbits. For groups
   with multiple non-identity operators, also pair atoms with each
   non-identity op (loop over the operator list instead of a single op)
   and union the `site_symmetry` arrays of any atom that's fixed by
   more than one op.

For a non-Cartesian-aligned op (e.g. inversion through a non-origin
point as in the diamond example), express the op's `t` vector in
**fractional** coordinates of the **simulation box** (which is
`N * lattice_a` per side for an N-cell setup). The diamond generator
sets `t = (1/(4N), 1/(4N), 1/(4N))` -- this is `2 * (1/(8N))`,
where `(1/(8N))` is the lamda position of the inversion center for an
NxNxN cell.

Validating your symmetry file
-----------------------------

The shipped JSON schema file tools/json/fix-symmetry-schema.json can
be used to validate the symmetry file structure and types.

In addition, `fix symmetry` performs four init-time checks that catch
the most common setup errors before the integrator starts:

1. **Box family** -- the box must match the declared `lattice` family
   (cubic, hexagonal, etc.).
2. **Tag existence** -- every atom tag listed in the orbit table must
   actually exist in the system.
3. **Image consistency** -- each declared image atom must sit at
   `R*s_asym + t (mod lattice)` to within the `tol` argument.
4. **Wyckoff consistency** -- each Wyckoff asym atom must satisfy the
   stabilizer constraints `(R_k - I)*s + t_k == 0 (mod lattice)` at
   init, also within `tol`.

If any check fails the fix prints a targeted error identifying the
offending orbit, op, and atom tag.

Common pitfalls
---------------

* **LJ atom overlap on dense lattices**. Tersoff and other materials-
  potential lattices are denser than typical LJ packing. If you adapt
  the diamond or SiC generator to LJ, set `lattice_a` so the
  nearest-neighbor distance is at least `~2^(1/6) * sigma`. The diamond
  generator default `a = 2.6` puts NN at the LJ minimum (NN = `a*sqrt(3)/4`).

* **Inversion that swaps atom types**. Zinc-blende structures (SiC,
  GaAs, ZnSe, ...) are *not* inversion-symmetric -- inversion swaps the
  two FCC sublattices. Use a different point-group element instead
  (e.g. a 2-fold rotation, a 4-bar rotation, or a non-diagonal mirror).
  The `generate_sic_c2.py` script demonstrates the C_2(z) choice.

* **Wyckoff sites on every atom**. Cubic structures like rocksalt have
  inversion centers on *every* atom site, so P-1 with inversion through
  origin makes all atoms Wyckoff with rank-3 constraints -- meaning
  every atom is pinned and the dynamics is trivial. Either shift to a
  non-symmetry-preserving lattice spacing, or pick a lower-symmetry
  subgroup that leaves more degrees of freedom free.

* **Off-origin inversion centers**. If your inversion center is not at
  the simulation-box origin (corner), the `t` vector of the op is
  non-zero and depends on the box size. The diamond example shows the
  required arithmetic. Express the center in fractional simulation-box
  coordinates and use `t = 2 * center`.

* **Atom numbering**. The generators write data files with explicit tags
  matching the basis-and-cell iteration order. If you switch from
  `read_data` to `create_atoms` with `lattice` + `region box`, the tag
  ordering can differ; regenerate the JSON with the new ordering. The
  generators use the same basis-cell iteration as LAMMPS's
  `create_atoms` for compatibility.

References
----------

* Cox, S. R. & White, A. D. *Symmetric Molecular Dynamics.* J. Chem.
  Theor. Comput. **18**, 4077-4088 (2022). DOI 10.1021/acs.jctc.2c00401
* Reference symd engine: https://github.com/whitead/symd
