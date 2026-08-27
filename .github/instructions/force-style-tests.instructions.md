---
applyTo: "unittest/**"
---

# LAMMPS Unit-Test Conventions (force-style YAML tests and friends)

Unit tests are CTest-based; build with `-D ENABLE_TESTING=on`, run with
`cd build && ctest -V [-R <pattern>]`.  Tests are organized by category under
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
