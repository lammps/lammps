# rRESPA Integration in Constraint/Force Fixes

Constraint fixes (`shake`, `rattle`, `ilves`) and forcing fixes that participate in
the inner-force update must support `run_style respa`.  The pattern is the same in all
of them; `fix_shake.cpp` is the canonical reference.  Steps to add respa support:

## If respa support is not feasible: refuse, never run silently wrong

Respa invokes only the `*_RESPA` fix hooks; masks without a respa variant
(`PRE_FORCE`, `PRE_REVERSE`, ...) are silently skipped, so a fix relying on them
does NOTHING under `run_style respa` (the ELECTRODE fixes silently froze their
electrode charges this way).  When implementing the hooks is out of scope, error
out in `init()`:

```cpp
if (utils::strmatch(update->integrate_style, "^respa"))
  error->all(FLERR, Error::NOLASTLINE, "Fix {} is not compatible with run_style respa", style);
```

Document the restriction in the fix's doc page, and if the fix has force-style YAML
tests, add it to the respa exclusion regex in `test_fix_timestep.cpp`.

## Header changes

- Forward-declare `class FixRespa;`.
- Add `void post_force_respa(int, int, int) override;` (and any other per-level hooks).
- Add private state: `int respa; int nlevels_respa; int *loop_respa; double *step_respa;
  class FixRespa *fix_respa; double dtf_inner;`.

## Source changes

- `#include "fix_respa.h"` and `#include "respa.h"`.
- Initialize the new members (especially `respa = 0`, `fix_respa = nullptr`) in the
  constructor.
- Add `mask |= POST_FORCE_RESPA;` in `setmask()`.
- In `init()`, detect respa with `utils::strmatch(update->integrate_style, "^respa")`;
  cache `Respa*` (for `nlevels`, `loop[]`, `step[]`) and `FixRespa*`
  (`modify->get_fix_by_style("^RESPA")` for `f_level[i][jlevel]`).
- In `setup()`, branch on `respa`: use `dtv = step_respa[0]`, set
  `dtf_inner = 0.5 * step_respa[0] * force->ftm2v` for the one-shot constraint
  projection, loop over levels calling `respa_ptr->copy_flevel_f(ilevel)`, your
  `post_force_respa(vflag, ilevel, loop_respa[ilevel]-1)`, then
  `respa_ptr->copy_f_flevel(ilevel)`.  Restore
  `dtf_inner = step_respa[0] * force->ftm2v` for normal stepping.
- In `reset_dt()`, branch on `integrate_style` so the inner-timestep variables stay
  consistent on a `timestep` reset.
- In your `unconstrained_update_respa(ilevel)`, set
  `dtfsq = dtf_inner * step_respa[ilevel]` and fold the inner-level contributions:
  `xshake = x + dtv*v + dtfsq/m*f + sum_{j<ilevel} 0.5*dt0*dt_j/m *
  fix_respa->f_level[i][j]`.

## Friend access to FixRespa::f_level

That array is `private` in `src/fix_respa.h`.  To read it from a new fix, add
`friend class FixYours;` alongside the existing `friend class FixShake;` /
`friend class FixRattle;` / `friend class FixIlves;` lines.

## Virial accumulation across respa levels (subtle, important)

Do NOT call `ev_init()` in `post_force_respa` -- it zeros the per-atom virial,
destroying the contribution accumulated at previously visited levels.  The correct
sequence (mirroring `FixShake::post_force_respa`):

1. Call `v_init(vflag)` *once*, at the innermost level on its final loop iteration.
2. Set `evflag = 1` only on each level's final loop iteration
   (`iloop == loop_respa[ilevel]-1`); set `evflag = 0` otherwise so intermediate
   iterations do not accumulate virial.

The framework relies on this exact discipline; getting it wrong shows up as stress
mismatches of order unity against the reference YAML even though forces and energies
match.

## Testing

The `unittest/force-styles/test_fix_timestep.cpp` driver exercises both verlet and
respa code paths for every YAML reference (respa with a `100 * epsilon` tolerance
multiplier); see `.github/instructions/force-style-tests.instructions.md`.
