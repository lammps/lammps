# unittest/granular — DEM verification tests

Verification tests for LAMMPS' granular (DEM) contact models. Each test drives a
small, analytically tractable system entirely from a YAML file and compares the
simulated result against a closed-form solution, and (once recorded) against a
stored per-atom reference trajectory.

## Layout

- `test_dem_01..15.cpp` — the fifteen test drivers. Each is a thin GoogleTest
  fixture (`newton_on` / `newton_off`) that hands off to the shared trajectory
  runner; all scenario-specific behavior lives in the YAML files.
- `test_dem_common.{cpp,h}` — builds the system from a YAML config, runs it in
  segments, and compares positions, velocities, torques, and angular velocities
  against the reference plus the selected analytic model.
- `test_analytic_models.{cpp,h}` — the closed-form models selected per test via
  the YAML `analytic_model` key.
- `test_config_reader.{cpp,h}`, `test_main.{cpp,h}`, `test_config.h` — YAML
  parsing and the test-runner entry point.
- `tests/` — the YAML test definitions, `demNN-<variant>-<dim>-<units>.yaml`.
- `dem_audit.py` — separate diagnostic script that compares results from all
  test drivers to analytic expectations. Intended for authoring/debugging
  new/existing tests. Not connected to CMake or CTest.

## Tests

Each numbered driver runs every `demNN-*.yaml` in `tests/` (registered
automatically by the CMake glob), so a driver may cover several variants of one
scenario across different contact models, dimensions, and parameters.

| # | scenario | analytic model |
|---|---|---|
| 01 | two-sphere head-on normal collision | `collision_restitution` |
| 02 | two-sphere elastic Hertzian normal impact (peak force) | `hertz_normal_impact` |
| 03 | sphere–wall elastic Hertzian normal impact (peak force) | `hertz_normal_impact` |
| 04 | oblique impact on a wall (gross-sliding friction) | `oblique_impact` |
| 05 | sphere sliding then rolling without slipping | `slip_cessation` |
| 06 | spinning sphere impact (rebound with friction) | `spin_impact` |
| 07 | rolling-resistance decay | `rolling_decay` |
| 08 | cohesive DMT pull-off force | `pulloff_dmt` |
| 09 | terminal velocity under fluid drag | `terminal_velocity_linear`, `terminal_velocity_schiller_naumann` |
| 10 | exact integration + static contact (free fall, stack) | `freefall`, `stack_energy` |
| 11 | contact with region walls | `wall_restitution` |
| 12 | superellipsoid collision | `momentum_conservation` |
| 13 | spinning sphere damped by twisting friction | `twist_decay`, `twist_decay_marshall` |
| 14 | granular heat conduction in static contact | `heat_equilibration` |
| 15 | two-sphere oblique impact (gross-sliding friction) | `oblique_impact_pair` |

Variants exercise the `hooke`, `hooke/history`, `hertz`, `hertz/material`,
`mindlin`, `mindlin/rescale`, and the `mindlin[_rescale]/force` models, the
`velocity`, `mass_velocity`, and `tsuji` damping variants, a 2D/LJ-units case,
and friction/angle
sweeps.  Variants named `demNN-legacy-*` (YAML tag `legacy`, also a CTest
label) exercise the classic `gran/hooke`, `gran/hooke/history`, and
`gran/hertz/history` pair styles and the classic `fix wall/gran` models.

## Building

Configure LAMMPS with testing and the GRANULAR package enabled, then build the
drivers:

    cmake -C ../cmake/presets/most.cmake -D ENABLE_TESTING=on -D PKG_GRANULAR=on ../cmake
    cmake --build . --target test_dem_01 test_dem_02 test_dem_03 test_dem_04 \
                            test_dem_05 test_dem_06 test_dem_07 test_dem_08 \
                            test_dem_09 test_dem_10 test_dem_11 test_dem_12 \
                            test_dem_13 test_dem_14 test_dem_15 -j

## Running

Run a single test file:

    ./test_dem_02 tests/dem02-hertzmaterial-twosphere-3d-si.yaml

Or run the whole suite through ctest from the build directory:

    ctest -R DEM --output-on-failure

Options:

- `-v` — verbose (show LAMMPS screen output).
- `-s` — print per-quantity error statistics.
- `-g <file>` — regenerate the reference into a new YAML file.
- `-u` — regenerate the reference in place.
- `-d <folder>` — set the folder for any external input files.

## YAML format

Each file specifies the system and the checks:

- `prerequisites`, `pre_commands`, `pair_style`, `pair_coeff`, `post_commands` —
  the LAMMPS commands that build the geometry, contact model, and integrator.
- `variables` — `${var}` substitutions used throughout the command blocks.
- `run_segments` — space-separated step counts; the trajectory is compared after
  each segment.
- `analytic_model`, `analytic_segment`, `analytic_tol` — which closed-form model
  to check, at which segment boundary, and to what relative tolerance.
- `run_pos` / `run_vel` / `run_torque` / `run_omega` — the recorded per-atom
  reference, produced by `-g` / `-u`.

## Regenerating reference data

After editing a test, regenerate its reference trajectory and confirm it passes:

    ./test_dem_0N tests/<file>.yaml -u
    ./test_dem_0N tests/<file>.yaml -s

The reference is recorded with `newton on`; tests are checked under both
`newton on` and `newton off`.
