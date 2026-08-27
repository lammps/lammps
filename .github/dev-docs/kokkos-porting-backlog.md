# KOKKOS Porting Backlog

**STATUS TABLES GO STALE.** Trust `ls src/KOKKOS/` and `git log` over this file.
Last verified: 2026-07-06.

History: pair Groups 1-4 (EXTRA-PAIR simple pairwise, Wolf/DSF/cut Coulomb, long-range
Coulomb; ~30 styles) and fix Groups A-C (trivial force loops, wall styles, NH-sphere
integrators, plus `compute temp/sphere`) were ported and merged onto `develop` in
May-June 2026. Nine fake (sync-to-host) fix ports were committed and deleted in the same
period ("remove fake KOKKOS ports" commit); see `git log` for details.

Coverage snapshot at last verification: pair coverage essentially complete (~101 styles;
all common ones done); kspace = plain `pppm` only; regions = `block` + `sphere` only.

## Open pair-style groups

- **Group 5 -- TIP4P** (high complexity: virtual M-site handling in the neighbor loop):
  `tip4p/cut`, `tip4p/long`, `lj/cut/tip4p/cut`, `lj/cut/tip4p/long`, plus kspace
  `pppm/tip4p`. **Ported on branch `kokkos-tip4p`** (pushed to the lammps/lammps remote):
  shared header-only base `PairTIP4PKokkos<DeviceType,PairCPUBase>` plus
  `PPPMTIP4PKokkos`. **PR not yet opened.** (Third-party TIP4P PRs #4971 and #4928 are
  closed.)
- **Group 6 -- granular/colloidal:** `gran/hooke` (no history arrays),
  `gran/hertz/history` (needs shear-history Kokkos views).
- **Group 7 -- MANYBODY three-body** (manual short-neighbor + three-body kernels; follow
  `pair_sw_kokkos` / `pair_tersoff_kokkos`): `sw/mod`, `nb3b/harmonic`, `tersoff/mod/c`.
- **Group 8 -- many-body / geometry-dependent** (custom kernels; the ASPHERE ones need
  per-atom quaternions on device): `atm` (triple loop), `local/density` (two-pass loop),
  `gayberne`, `resquared`, `line/lj`, `tri/lj`.
- **Group 9 -- H-bond, dipole, per-atom sphere:** `hbond/dreiding/lj`,
  `hbond/dreiding/morse`, `hbond/dreiding/lj/angleoffset`,
  `hbond/dreiding/morse/angleoffset` (angle-dependent triple loops), `lj/sf/dipole/sf`,
  `lj/expand/sphere`.

## Fix/compute work in open PR #4989 (branch `more-kokkos-porting`, not yet merged)

- Group D (per-atom stored state, KokkosBase): `fix spring/rg`, `ti/spring`,
  `addtorque/group`.
- Group E (thermostat/barostat coupling): `fix gjf`, `press/berendsen`,
  `press/langevin`, `temp/csvr`, `temp/csld`.
- Computes: `temp/partial`, `temp/profile`, `ke/atom`.

## Group F triage

- Real porting targets: `fix heat` (group-KE `parallel_reduce` + velocity-rescale
  kernel), `fix indent` (per-atom force kernel; host-evaluate only the few variable
  scalars per step and broadcast them into the kernel).
- LEAVE UNPORTED (full-ports-only policy; no device work to offer): `fix move`
  (`xoriginal` is a plain `double **` CPU array; a real port means converting it to a
  DualView + KokkosBase first), `fix restrain` (`atom->map()` ghost lookups host-side),
  `fix deform/pressure` (box-deformation logic is host-side scalar work).

## Priority order (metric: device residency; NPT/barostatting first)

1. Land the open work: PR #4989, and open a PR for branch `kokkos-tip4p` (NPT of TIP4P
   water is a flagship GPU use case; the per-step host pair+kspace fallback -- not the
   barostat -- is the real bottleneck).
2. `fix heat/kk` and `fix indent/kk` (the Group F real targets).
3. `compute temp/region/kk` -- gated on region coverage (only `block` and `sphere`
   regions are ported; port the needed `region_*_kokkos` first).
4. `fix box/relax/kk` (minimization-time barostat; lower frequency).
5. Pair Groups 6-9, case by case.

Deferred (not on any branch):

- `ewald/kk`: Ewald is O(N^3/2); GPU pays off for large systems, where `pppm/kk` is what
  is actually used. Ewald's real role is a CPU validation reference for small test
  cases. The kspace GPU target is the PPPM family (e.g. `pppm/disp/kk`).
- `compute stress/atom`, `pe/atom`, `centroid/stress/atom`: need framework work first
  (see "per-atom virial/energy aggregation" in the porting guide); a host-aggregation
  wrapper would be a fake port.

## Do NOT port

- **`fix npt/kk` and `nph/kk` are COMPLETE full ports.** `fix_npt_kokkos` /
  `fix_nph_kokkos` create a `temp/kk` temperature and a `pressure` compute that *uses*
  it; every integrator hook (`nh_v_press`, `nve_v`, `nve_x`, `nh_v_temp`) and `remap()`
  (via `DomainKokkos::x2lamda`/`lamda2x`) runs as a device kernel, including the
  triclinic path. Nothing remains to do on the Nose-Hoover barostat itself.
- **Do NOT port `compute pressure`.** It has `datamask_read = datamask_modify =
  EMPTY_MASK`: no per-atom loop at all -- the scalar/vector is the temperature scalar
  (already from `temp/kk`) plus the global `virial[6]` arrays (already reduced on-device
  by the force styles). The `if (!pressure->kokkosable) atomKK->sync(...)` guard in
  `fix_nh_kokkos` therefore syncs **nothing**. A `pressure/kk` wrapper would eliminate a
  no-op and add zero GPU work.
- Any fix whose per-atom work cannot run on-device: `move`, `restrain`,
  `deform/pressure` (see the Group F triage above), and anything coupling to an external
  library (`colvars`). Leaving them unported is the correct outcome.
- `fix ilves/omp` -- per-cluster OpenMP threading was benchmarked and rejected (net
  slowdown); the same reasoning argues against a KOKKOS thread-parallel variant.
