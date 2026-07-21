# unittest/granular replacement — mapping notes

Only the 8 tests derived from your scripts are here, renumbered 01-08
(sequentially, in original scenario order). Everything else has been
removed — no leftover placeholders for unimplemented categories, no
unused test_dem_0N.cpp files.

## What's unchanged (kept from the existing directory)

Generic test *harness* files — not test content, and match LAMMPS' own
documented design for `unittest/granular`:

- `CMakeLists.txt` — trimmed to register exactly 8 executables/test sets
- `test_config.h`, `test_config_reader.{h,cpp}`
- `test_dem_common.{h,cpp}`, `test_main.{h,cpp}`
- `test_analytic_models.{h,cpp}` — already implements every closed-form
  model these 8 tests need (`slip_cessation`, `oblique_impact`,
  `rolling_decay`, `pulloff_dmt`, `collision_restitution`, `spin_impact`)

## Script → test mapping (new numbering)

| # | Your script | New file | analytic_model |
|---|---|---|---|
| 01 | `sphere_sliding_rolling.in` | `dem01-hertzmindlin-cof-3d-si.yaml` | `slip_cessation` |
| 02 | `oblique_impact.in` | `dem02-oblique-cof-3d-si.yaml` | `oblique_impact` |
| 03 | `sphere_pure_rolling.in` | `dem03-rolling-cof-3d-si.yaml` | `rolling_decay` |
| 04 | `in.cohesive_dmt` | `dem04-dmt-twosphere-cof-3d-si.yaml` | `pulloff_dmt` |
| 05 | `angular_sphere_impact.in` | `dem05-spin-headon-cof-3d-si.yaml` | `collision_restitution` |
| 06 | `in.angular_bidisperse` | `dem06-spin-unequal-cof-3d-si.yaml` | `spin_impact` |
| 07 | `normal_contact_cof.in` | `dem07-wall-normal-cof-3d-si.yaml` | `hertz_normal_impact` |
| 08 | `in.normal_simple` | `dem08-twosphere-normal-3d-si.yaml` | `hertz_normal_impact` |

Two of these (02, 06) have their geometry relabeled — not their physics
changed — to match axis conventions `test_analytic_models.cpp` hard-codes
(e.g. `oblique_impact` reads `v[i][0]`/`v[i][2]` specifically). Every
material property, density, diameter, velocity, friction coefficient, and
restitution coefficient is transcribed unchanged from your original
scripts. See the previous round's notes below for exact per-file details.

## What you still need to do — this is important

**I could not run LAMMPS in this environment.** Each of the 8 YAML files
has a real `prerequisites`/`pre_commands`/`variables`/`pair_style`/
`pair_coeff`/`post_commands` block built from your script, plus a
placeholder `run_segments` guess, but the `run_pos`/`run_vel`/
`run_torque`/`run_omega` reference blocks are marked `# TODO: regenerate`.

Generate the real reference data by building LAMMPS with `GRANULAR` +
testing enabled and running, e.g.:

```
test_dem_01 tests/dem01-hertzmindlin-cof-3d-si.yaml -g new.yaml
```

or update in place with `-u`, then check with `-s` (per-quantity error
stats) that `analytic_segment` lands where the model expects (e.g. peak
compression for `hertz_normal_impact`, post-rebound free flight for
`oblique_impact`/`spin_impact`) and retune `run_segments` if not.

Per-file notes carried over from before:
- `dem01` (from `sphere_sliding_rolling.in`) and `dem03` (from
  `sphere_pure_rolling.in`): source scripts used gravity magnitude `1.0`
  (not standard 9.81) — kept as literally scripted, verify that's intended.
- `dem01`: source script had a 10,000-step settle phase before applying
  sliding velocity; this version places the sphere directly at contact
  instead, so the early transient will differ slightly from your original.
