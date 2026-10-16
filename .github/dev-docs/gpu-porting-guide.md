# GPU Package Porting Guide (pair styles)

How to add a `/gpu` variant of a pair style to the GPU package.  Written from the
2026-10 batch of 18 charge-family ports (`coul/wolf`, `lj/charmmfsw/coul/long`, the
FEP soft-core and CORESHELL `/cs` styles, ...); clone from those, they use the
current API.  The GPU package has two halves: the device library `lib/gpu/` (Geryon
`UCL_*` containers, one kernel source compiled for CUDA, HIP, and OpenCL) and thin
LAMMPS wrappers in `src/GPU/`.

## Files of one port (6 files + 3 edits)

| file | role |
|---|---|
| `lib/gpu/lal_<name>.h`, `.cpp` | device class `<Name><numtyp,acctyp>` derived from a base class; uploads per-type parameters, launches the kernels in `loop()` |
| `lib/gpu/lal_<name>.cu` | kernels `k_<name>` (generic) and `k_<name>_fast` (per-type parameters in shared memory) |
| `lib/gpu/lal_<name>_ext.cpp` | C++ entry points `<prefix>_gpu_init/_reinit/_clear/_compute_n/_compute/_bytes` around one static device object (`<PREFIX>MF`) |
| `src/GPU/pair_<style>_gpu.h`, `.cpp` | `Pair<Style>GPU : public Pair<Style>`; calls the entry points |
| `lib/gpu/include/lammps_gpu.h` | declarations of the entry points in namespace `LAMMPS_GPU` |
| `src/GPU/Install.sh` | two `action` lines per style |
| `doc/src/` | `(g)` letter in `Commands_pair.rst`, `.. index:: pair_style <style>/gpu`, `*<style>/gpu*` in the "Accelerator Variants" line |

- The kernel symbol and generated header name are the `.cu` file name minus `lal_`:
  `lal_coul_wolf.cu` -> `coul_wolf_cl.h` / `coul_wolf_cubin.h`, string `coul_wolf`.
  The `.cpp` includes the matching header under `USE_OPENCL` / `USE_CUDART` / else.
- CMake globs `lib/gpu/*.cu` and `*.cpp` (`cmake/Modules/Packages/GPU.cmake`):
  reconfigure (`cmake -S cmake -B <build>`) after adding files.  A kernel that includes
  an extra header (`lal_*_extra.h`) must be removed from the OpenCL glob list and get
  its own `GenerateOpenCLHeader` line there, like `tersoff` and `gayberne`.
- No `src/Purge.list` entry: that file lists REMOVED or renamed files only.
- `Install.sh`: `action pair_<style>_gpu.cpp pair_<style>.cpp` (and the same for `.h`),
  where the second file is the CPU parent, so the variant is installed only when the
  parent's package is.  Keep the block alphabetical.
- `lammps_gpu.h`: copy the block of the template style, then format only the new lines
  with `clang-format -i --lines=<first>:<last> lib/gpu/include/lammps_gpu.h` (uses
  `lib/gpu/include/.clang-format`; do not reformat the whole file).  A declaration
  that does not match the `_ext.cpp` definition becomes a LINK error, not a compile
  error, because the definition then declares a new overload.

## Choosing the template

Clone all six files of the closest existing style (a sed/python rename of class name,
`MF` object, function prefix, kernel names, header guard, include names), then change
the physics.  Base class: `BaseAtomic` (no charges), `BaseCharge` (needs `q`; use it
for anything with Coulomb), `BaseDipole`, `BaseDPD`, `BaseThree` (sw/tersoff
family), `BaseEllipsoid`, `BaseSPH`, `BaseAmoeba`.

| new style needs | clone from |
|---|---|
| LJ-like with offset, per-type params | `lal_lj96` (NOT `lal_lj`: its `ONETYPE` specialization adds complexity a clone does not need) |
| cut Coulomb / LJ + cut Coulomb | `lal_coul`, `lal_lj_coul`, `lal_lj_class2_coul_cut` |
| Ewald real space (erfc polynomial) | `lal_coul_long`, `lal_lj_coul_long`, `lal_buck_coul_long` |
| Wolf / DSF | `lal_coul_wolf`, `lal_coul_dsf`, `lal_lj_cut_coul_wolf`, `lal_born_coul_dsf` |
| CHARMM switching, force switching | `lal_charmm_long`, `lal_charmmfsw_long`, `lal_charmmfsw_fsh` |
| core/shell `/cs` | `lal_coul_long_cs`, `lal_born_coul_long_cs`, `lal_coul_wolf_cs` |
| FEP soft core (`lambda` per type pair) | `lal_lj_cut_soft`, `lal_lj_coul_soft`, `lal_lj_coul_long_soft`, `lal_coul_cut_soft` |
| global parameters only | `lal_coul_slater_cut` (scalar kernel arguments) |

Long-range Coulomb on the GPU is ALWAYS the analytic erfc polynomial: there are no
Coulomb tables, so the `*_table*` YAML tests of a new long-range style get `gpu` in
`skip_tests`.  KSpace stays on the host (or uses `pppm/gpu`).

## Device class (`lal_<name>.h/.cpp`)

- `init()` calls `init_atomic(nlocal, nall, max_nbors, maxspecial, cell_size, screen,
  <name>, "k_<name>")` (same name in `BaseCharge`), then decides `shared_types`: when
  `ntypes <= device->max_shared_types()` and the block is big enough, the per-type
  arrays are padded to `max_shared_types` and the `_fast` kernel runs.  The YAML decks
  (4-5 types) therefore test the FAST kernel; the generic kernel runs only with many
  types or a small block size.
- Per-type parameters: `UCL_D_Vec<numtyp4>` (or `numtyp2`, `numtyp`) filled with
  `atom->type_pack4/2/1(ntypes, lj_types, dev, host_write, a, b, c, d)` from the host
  `double **` arrays; global scalars become kernel arguments (`numtyp` members); a
  handful of derived constants can go into a small `UCL_D_Vec<numtyp>` (see the `fsw`
  array of `lal_charmmfsw_long`).  Special-bond factors: `sp_lj` with 4 LJ (+4 Coulomb)
  entries.
- `reinit()` re-packs the per-type arrays with the stored `_lj_types`.  Provide it for
  every style whose parameters `fix adapt` can change (`extract()` in the CPU style);
  otherwise the wrapper must disable reinit (see the host wrapper section).
- `loop()` launches `k_pair_sel` (fast) or `k_pair` (generic) with `run(&a, &b, ...)`.
  The argument list must match the kernel signature EXACTLY in order and count, for
  BOTH kernels (the generic one has the extra `lj_types` argument).  No back end checks
  this at compile time; a mismatch shows up as wrong numbers or a launch failure.

## Kernels (`lal_<name>.cu`)

- Types: `numtyp` = float (single, mixed) or double; `acctyp` = accumulator, double in
  mixed and double precision.  Write literals as `(numtyp)1.0`; use `ucl_recip`,
  `ucl_rsqrt`, `ucl_sqrt`, `ucl_exp`, `ucl_powr` so all three back ends compile.  Small
  helpers: `ucl_inline numtyp f(...)` (OpenCL C has no references or templates).
- Skeleton (copy it): `atom_info`, `local_allocate_store_charge()` (or `_pair`),
  `nbor_info`, the neighbor loop with `sbmask(j)` / `j &= NEIGHMASK`, then
  `store_answers_q(...)` (or `store_answers`).  Neighbor lists are FULL and only local
  atoms receive forces, so results do not depend on the newton setting.
- Energies: every pair is visited twice and `store_answers*` multiplies by 0.5.  Per-atom
  terms that must be counted ONCE (Wolf/DSF self energy) are added as
  `e_coul += (acctyp)2.0*e_self` with `e_self` divided by `t_per_atom`, as in
  `lal_coul_dsf.cu`.
- Special bonds: cut Coulomb uses `factor_coul = sp_lj[sbmask(j)+4]`; erfc-based
  (long, DSF, `/cs`) kernels use `factor_coul = (numtyp)1.0 - sp_lj[...]` and subtract
  `factor_coul*prefactor`, matching the CPU `(1.0-factor_coul)*prefactor` term.  Check
  which convention the template uses before changing the formula.
- Constants: `EWALD_F`, `EWALD_P`, `A1`..`A5`, `MY_PIS` come from `lal_preprocessor.h`;
  the `/cs` kernels define their own `CS_EWALD_P`, `B0`..`B5`, `CS_EWALD_F`, `EPSILON`
  and do `rsq += EPSILON` inside the cutoff, with `acctyp` for the erfc terms and
  `erf()` for excluded pairs.
- Fast kernel: copy into `__local` arrays EVERY per-type array that the FORCE uses.
  The template pattern `if (EVFLAG && eflag) lj3[tid]=lj3_in[tid];` is only valid when
  `lj3` feeds the energy alone; soft-core styles use `lj3` in the force, and the
  eflag-only copy read uninitialized local memory on steps without energy tally
  (fixed in `lal_lj_coul_soft.cu` / `lal_lj_coul_long_soft.cu` in 2026; the YAML gpu
  tests of those styles had been skipped as "unstable").
- Global LJ cutoffs (CHARMM): pass `cut_ljsq`, `cut_lj_innersq`, `denom_lj` as kernel
  arguments; do not abuse the per-type `cutsq` slot.

## Entry points (`lal_<name>_ext.cpp`)

- `init_message(screen, "<style>", first_gpu, last_gpu)` prints "Using acceleration for
  <style>": use the exact style name without `/gpu` (a cloned file keeps the template's
  name otherwise; 22 styles printed a wrong one until 2026-10).
- Keep the init sequence of the template: world rank 0 initializes and compiles first,
  then the other ranks sharing a device, with `world_barrier()`/`gpu_barrier()`.  The
  same rank loop applies to `_reinit`.

## Host wrapper (`src/GPU/pair_<style>_gpu.{h,cpp}`)

- Constructor: `respa_enable = 0; suffix_flag |= Suffix::GPU;
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);`.  Set `reinitflag = 0` when the
  wrapper does NOT implement `reinit()`: `fix adapt` then stops with an error instead
  of changing only the host copy of the parameters.
- `compute()`: two paths.  `gpu_mode != GPU_FORCE` -> `_compute_n` builds the
  neighbor list on the device (needs sublo/subhi, tag, special); `GPU_FORCE` -> the
  host full list (`list->ilist` ...) goes to `_compute`.  Afterwards
  `if (atom->molecular != Atom::ATOMIC && neighbor->ago == 0) neighbor->build_topology();`.
- `init_style()`: repeat the `cutsq` loop with `init_one(i,j)` (Pair::init calls
  `init_one` only AFTER `init_style`), compute every derived constant the CPU style
  computes in `init_style`/`init_one` (switching denominators, soft-core factors,
  force-switch constants), call `<prefix>_gpu_init(...)`, `GPU_EXTRA::check_flag()`,
  and request a full neighbor list only in `GPU_FORCE` mode.
- `reinit()`: `Pair::reinit();` (recomputes via `init_one`) then `<prefix>_gpu_reinit()`
  with the updated arrays.
- Hybrid pair styles need host neighbor lists (`-pk gpu N neigh no`); `fix gpu` stops
  with an error otherwise.
- Do not copy a CPU bug into the kernel: fix the CPU style (and its `/omp`, `/kk`
  variants) in a separate commit.  Example: `morse/soft` `compute()` never subtracted
  the `pair_modify shift` offset.

## Testing

- Build double precision first (`-D GPU_PREC=double`): it gives tight comparisons.
  Then mixed (the default).  Since the CMake rework (#5198) CUDA architectures are
  numbers: `-D GPU_ARCH=70` (the old `sm_70` is converted with a deprecation warning).
- Compile-only check without a GPU run: a HIP build compiles all `.cu` files with
  `cmake --build <build> --target gpu`; OpenCL kernels compile only at run time, so an
  OpenCL run is required before calling a kernel portable.
- Run new kernels on a dedicated test machine when possible: a faulting kernel can
  hang a GPU that also drives the desktop.
- YAML (`unittest/force-styles/tests/`): the `gpu` case reuses the CPU reference with
  `-sf gpu`, newton off, device neighbor lists (`neigh no` for `hybrid*` styles or the
  `host_neigh` tag), and epsilon x7.5 (double), x5e8 (mixed), x1e10 (single).  Tags
  `gpu_no_single`/`gpu_no_mixed` skip a precision.  A missing `/gpu` variant or a skip
  makes the case pass silently: confirm it RAN (`ctest -V`, `lmp -h`).  When a mixed
  precision failure is a single near-zero force component (relative error of the
  largest component ~1e-6), relax the YAML epsilon instead of changing the kernel.
- MD comparison: 20 steps of each YAML system plain vs `-sf gpu` with device and host
  neighbor lists (host only for hybrid styles), newton on, and 2 MPI ranks sharing one
  GPU; compare thermo (energies, pressure tensor, summed |f|) and per-step sums of
  `pe/atom` and `stress/atom`.  Expect 1e-12 or better in double precision, about 1e-6
  to 1e-5 in mixed.  For styles with `reinit()`, add a deck with `fix adapt` ramping a
  parameter every step: it is the only test of the reinit path.  Also run a deck where
  the mixing rule holds (no explicit off-diagonal `pair_coeff`) if the fast kernel
  mixes on the fly.  The decks and scripts in `kokkos-tier1-verification/` show the
  pattern for such a comparison driver.
