---
applyTo: "unittest/**"
---

# LAMMPS Unit-Test Conventions (force-style YAML tests and friends)

Unit tests are CTest-based; build with `-D ENABLE_TESTING=on`, run with
`ctest --test-dir build -V [-R <pattern>]`.  Tests are organized by category under
`unittest/` (`force-styles/`, `commands/`, `formats/`, `c-library/`, `fortran/`,
`python/`, `utils/`, `granular/` -- the latter has its own instructions file).

## YAML-driven force-style tests (`unittest/force-styles/`)

- Each test driver (`test_pair_style`, `test_bond_style`, `test_angle_style`,
  `test_fix_timestep`, ...) loads a `.yaml` reference file and compares thermo,
  forces, energies, and stresses against it with `epsilon`-scaled tolerances.
- The drivers automatically exercise EVERY accelerator-suffix variant of a style
  (`/omp`, `/intel`, `/gpu`, `/kk`) that is compiled into the test executable, from
  the single base-style YAML reference.  Adding an accelerated variant needs NO new
  YAML file -- just run the existing reference with that package enabled.
- Regenerate or update reference data with the driver's command-line flags
  (`-g <file>` generate, `-u` update in place, `-s` print per-quantity error
  statistics for tuning `epsilon`).  Prefer `-u` so the file history stays clean.
  Regeneration can reset the `tags:` line -- re-check it after every `-u`.
- A YAML with a missing prerequisite or `input_coeffs` entry SKIPS silently while
  ctest still reports "Passed".  After adding or editing a YAML, confirm from the
  gtest output that its cases actually executed.
- The restart leg re-applies `post_commands` after `read_restart`, and restart files
  do not store `neigh_modify` settings -- put those into `post_commands` when a test
  depends on them.  A fix whose random-number state is not written to restart files
  cannot pass the restart leg; store the state (`RanMars::get_state()`/`set_state()`)
  rather than skipping the test.

## Tolerances and portable reference data

- `epsilon` is the relative tolerance of the plain runs; accelerator sub-tests scale
  it (e.g. in `test_pair_style`: GPU x7.5 double, x5e8 mixed, x1e10 single; KOKKOS
  x5, with a further x2e9 for mixed and x1e10 for single precision builds).
- Raising `epsilon` to about 5e-13 is acceptable for analytical kernels once `-s`
  shows the residual is precision noise (compare with the error profile of a
  sibling style on the same input).  Beyond 1e-12 for an analytical style needs
  maintainer approval -- it also loosens the CPU/OPENMP/KOKKOS comparisons and hides
  real bugs.  Spline- and table-based styles are legitimately noisier.
- Write references that do not depend on the last bits: define `fix enforce2d` after
  every fix that adds forces or torques, use `velocity ... loop all` (not
  `loop geom`), and avoid reference quantities that are pure roundoff.  Before
  loosening a tolerance or tagging a test `unstable` for an ARM64- or macOS-only
  failure, follow the platform triage in
  `.github/dev-docs/testing-and-verification.md`.

## Torque coverage

- Per-atom torque trajectories are recorded and compared ONLY by the
  `test_fix_timestep` driver (`run_torque` blocks).  To lock in the torque behavior
  of a pair style (e.g. dipole styles), add a fix-timestep fixture with
  `fix ... nve/sphere update dipole` and that pair style in `post_commands`
  (precedent: `fix-timestep-nve_sphere_dipole_ljlong.yaml`, which pins the
  LJ-only cutoff shell where a torque bug once hid).
- The `ellipsoid` entry on a `tags:` line makes `test_pair_style` ALSO assert
  `pair->single()` extra output (`svector` forces+torques, `single_extra >= 6`).
  Only tag styles that implement that interface (gayberne/resquared family) --
  never dipole styles.

## rRESPA coverage in fix tests

- `test_fix_timestep` exercises BOTH the verlet and respa code paths for every YAML
  reference; the respa path automatically applies a `100 * epsilon` tolerance
  multiplier.
- If a fix genuinely cannot support `run_style respa`, add it to the exclusion regex
  in `test_fix_timestep.cpp` -- but the strongly preferred fix is to make the fix
  respa-compatible; see `.github/dev-docs/respa-integration.md` for the standard
  pattern (including the subtle virial-accumulation discipline).
- Respa-related stress mismatches of order unity against the reference YAML, with
  matching forces and energies, almost always mean `ev_init()` was called in
  `post_force_respa()` (it must not be; see the respa guide).

## Adding tests for new styles

- Copy the YAML of a closely related style, adjust the input setup, regenerate the
  reference data with `-g`/`-u`, then re-run CMake so new files register with CTest.
- Verify new force styles against numerical differentiation (`fix numdiff`) where
  possible; see `.github/dev-docs/testing-and-verification.md`.
- The generators add a `generated` entry to the `tags:` line of newly written
  files on purpose: it marks reference data that has not been reviewed yet and
  makes such files easy to find with grep.  Remove the tag as the LAST step,
  after the reference data has been reviewed and validated -- never leave it in
  a file you commit as finished work.
