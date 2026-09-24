---
applyTo: "unittest/granular/**,src/GRANULAR/**"
---

# Granular / DEM Test Suite (`unittest/granular/`)

YAML-driven regression + analytic tests for DEM/granular models.  Eleven programs
`test_dem_01`..`test_dem_11` (43 tests): DEM-01..06 mirror the MFiX-DEM VVUQ cases;
DEM-07..11 add rolling resistance, cohesion (JKR/DMT pull-off), two-sphere and
spinning-sphere collisions, bulk angle-of-repose/rigid clumps, and elastic Hertzian
normal impact (benchmark sets from Mohajeri 2024, Saomoto 2023, and Chung & Ooi 2011).
Documented in `doc/src/Developer_unittest.rst` ("Tests for granular (DEM) models").

## Architecture

All driver logic lives in `test_dem_common.cpp` (compiled into the `granular_tests`
library); each `test_dem_0N.cpp` holds only the `newton_on`/`newton_off` GoogleTest
fixtures.  The whole system is built *from the YAML* (not a fixed input template): a
`variables` block becomes `${var}` index variables, `pre_commands` build the geometry,
`pair_style`/`pair_coeff` set the contact model, `post_commands` add fixes, and
`run_segments` advances the trajectory capturing per-atom `pos/vel/torque/omega/angmom`
after each segment.  Closed-form checks are opt-in via `analytic_enable/_model/_tol/
_segment` and implemented in `test_analytic_models.cpp`.  The harness shares
`yaml_reader.h`, `yaml_writer.{h,cpp}`, and `error_stats.{h,cpp}` with
`unittest/force-styles/` via `../` includes; `test_config*` and `test_main*` are
granular-specific forks.  For chaotic bulk tests (e.g. angle of repose) whose per-atom
trajectory is not bit-reproducible across newton on/off or platforms, set
`analytic_only: yes`: the generator records no per-atom blocks and only the analytic
observable is checked (DEM-10 pairs this with a short deterministic settling regression
plus a rigid-clump case so the bit-for-bit path stays covered).

## Adding a YAML reference

Copy a similar `dem0N-*.yaml`, edit variables/commands, then regenerate the reference
data **in place** with `TEST_ARGS=-u ctest -R DEM0N:variant` (or `test_dem_0N file.yaml
-u`).  Re-run CMake so the new test registers; the driver `-s` flag prints per-quantity
error statistics for tuning `epsilon`.

## Adding a test program

Thin-copy `test_dem_0N.cpp` (only the suite name changes), add an `add_executable`/
`register_dem_tests` pair to `CMakeLists.txt`, add `dem0N-*.yaml` files; if needed add
a named model to `test_analytic_models.cpp` (read params from `variables`, read masses
and radii from the live `atom` arrays, assert relative error with `EXPECT_LE`).

## Gotchas / lessons

- **Never regenerate to a `dem0N-*.yaml` sibling** (e.g. `-g new.yaml` in `tests/`):
  the `CONFIGURE_DEPENDS` glob registers it as a phantom/stale test until the next
  reconfigure.  Use `-u` (in place); idempotent apart from the date line.
- Granular/atomic systems have **no atom map** by default -> the reference generator
  iterates local atoms by `tag[i]`, never `atom->map(tag)` (segfaults).  `EXPECT_*`
  index by tag.
- Per-quantity flag guards (`torque_flag`/`omega_flag`/`angmom_flag`) auto-select
  fields: spheres carry `omega` (no `angmom`); ellipsoid/superellipsoid carry `angmom`
  (no `omega`).
- The libyaml reader accepts only scalar values (incl. block literals) -- all
  multi-value fields are newline/space-split block scalars; per-atom blocks use
  `segment tag x y z` rows.
- Analytic tolerances are loose (1e-3..5e-2): soft-sphere DEM only *approaches* the
  hard-sphere/instantaneous-contact ideal.  Use `coeff_restitution` when the analytic
  model needs a known restitution `e`.  Free-fall is exact (velocity-Verlet integrates
  constant acceleration exactly, ~1e-15).
- Soft contacts corrupt closed forms in two more ways: rigid-rolling kinematics
  (`omega = v/(r - delta)`) needs `delta/r << tol`, and in two-sphere oblique impacts
  the line of centers rotates during the finite contact, throwing fixed-normal impulse
  forms off by ~5% at `emod` 7e10.  Stiffen the contact (e.g. `emod` 7e11 / `kn` 3e8
  brings it to ~0.1%) instead of loosening tolerances.
- The harness requires `pair_style`/`pair_coeff` keys even for wall-only tests; a
  missing key produces a confusing error from the default `pair_style zero` fallback.
- Physics setup: start a particle resting on a wall at the gravity-equilibrium overlap
  (`z0 = r - mg/kn`) so the normal force is steady at `mg` (DEM-04).  The gross-sliding
  oblique-impact closed form needs a grazing angle (`vx_in > (7/2)mu(1+e)vz_in`)
  (DEM-05).  The modular `pair_style granular` tangential `linear_history` needs
  nonzero tangential damping or it rings undamped around the no-slip state;
  `gran/hooke/history` (dampflag 1) settles on its own.
- The historic `fix wall/gran hooke/history` energy-injection-at-grazing bug is fixed
  (develop commit `8285aa04b2`, `contact_radius_flag = 0` for the classic tangential
  model in `gran_sub_mod_tangential.cpp`).  Grazing/gross-sliding oblique tests run
  cleanly to ~76 deg (DEM-05); the 76-deg `dem05-hookehistory-grazing`
  `energy_dissipation` test locks the fix in.  Do NOT reintroduce the old "keep
  oblique tests <=45 deg" restriction.
- DEM-06 exercises the production fix `viscous/nonlinear` (Schiller-Naumann drag) in
  `src/GRANULAR/`; it reduces to Stokes drag (`6 pi mu r v`) at low Reynolds number.
- DEM-10/11 setup: `fix pour` needs the insertion region wider than one diameter AND
  `volfrac*vol/vol_one >= 1` (else "insertion count per timestep is 0"); a perfect
  lattice is mechanically stable and won't form a repose angle, so a real heap needs a
  seeded pour.  Rigid clumps use `atom_style hybrid sphere molecular` + `fix
  rigid/small molecule` and need `thermo_style custom step atoms` to dodge a spurious
  negative-DOF temperature error.  The Hertz peak-impact check (`hertz_normal_impact`)
  asserts the convention-independent energy balance
  `1/2 mu_red V^2 = 2/5 P_max alpha_max` (the 2/5 is the delta^1.5 signature); time a
  `run_segments` boundary to land at peak compression (relative normal velocity ~0).
