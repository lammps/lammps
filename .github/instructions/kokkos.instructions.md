---
applyTo: "src/KOKKOS/**"
---

# KOKKOS Package Rules (always-on)

These rules apply to every change under `src/KOKKOS/` and are the canonical single copy:
review code against them directly. How-to templates and lessons live in the deep-dive docs
listed at the end; build and test commands live in `.github/copilot-instructions.md`.

## Key rules

1. **`if (copymode) return;` must be the FIRST line of every base-class destructor** that a
   `_kokkos` style inherits from. A KOKKOS style object is copied by value for device
   execution (`copymode = 1`); without the guard, the copy's destructor frees memory owned
   by the original (double free / GPU crash). Check with
   `grep -n copymode src/<PKG>/<base>.cpp`.

2. **Declare `virtual void allocate()` in the base-class header.** The KOKKOS subclass
   overrides `allocate()` to replace plain arrays with dual views; without `virtual`, the
   base version is silently called and the dual-view setup never happens (wrong results or
   crashes).

3. **Per-atom views assigned from `atomKK->k_<field>.view<DeviceType>()` must be declared
   `typename AT::t_<type>`, never `DAT::t_<type>`.** `AT` is `ArrayTypes<DeviceType>` and
   adapts to the template parameter; `DAT` is hardcoded to `ArrayTypes<LMPDeviceType>`.
   GPU builds instantiate the style for BOTH `LMPDeviceType` and `LMPHostType` (the
   `kk/host` variant); a `DAT::` per-atom view then hits a Kokkos `ViewMapping` static
   assertion at compile time. On CPU-only builds both types share one memory space, so the
   bug is silent there and only surfaces when someone builds for GPU. `DAT::` is only safe
   for dual views the style allocates itself and never assigns from `atomKK` (e.g.
   `k_eatom`, `k_vatom`, `k_cutsq`, `k_cut_ljsq`, `k_cut_coulsq`). Audit with
   `grep "DAT::" src/KOKKOS/<style>_kokkos.h`: every hit must be a self-allocated dual view.

4. **Use `Kokkos::` math functions in device kernels** (`Kokkos::sqrt`, `Kokkos::exp`,
   `Kokkos::pow`, `Kokkos::log`, ...), never `std::` -- `std::pow` and friends fail to
   compile for GPU targets. `powint()` from `math_special.h` is host-only; replace it with
   explicit repeated multiplication in device code.

5. **Register new files in `src/KOKKOS/Install.sh`** with `action` lines: omit the second
   argument when the base class is in core `src/`; pass both filenames when the base is in
   an optional package:

   ```sh
   action fix_lineforce_kokkos.cpp             # base in src/
   action fix_drag_kokkos.cpp fix_drag.cpp     # base in EXTRA-FIX
   ```

## Datamask and sync/modify discipline

- In the constructor set `kokkosable = 1;`, `atomKK = (AtomKokkos *) atom;`,
  `execution_space = ExecutionSpaceFromDevice<DeviceType>::space;`, and declare
  `datamask_read` / `datamask_modify` with the bitmask constants from `atom.h` (`X_MASK`,
  `V_MASK`, `F_MASK`, `MASK_MASK`, ...): exactly the per-atom arrays the style reads and
  writes. Override `get_compute_flag()` if needed.
- Before a kernel reads per-atom data: `atomKK->sync<DeviceType>(datamask_read)` (or
  `atomKK->sync(execution_space, datamask_read)`), then assign the view members. After a
  kernel writes: `atomKK->modified<DeviceType>(datamask_modify)`.
- Fixes that keep per-atom state across atom migration (they implement `grow_arrays`,
  `copy_arrays`, `pack_exchange`, `unpack_exchange`) must convert those arrays to
  `Kokkos::DualView` members and also inherit from `KokkosBase`, overriding
  `pack_exchange_kokkos` / `unpack_exchange_kokkos` (pattern: `fix_spring_self_kokkos`).

## FULL PORTS ONLY (project policy)

We accept **only full ports**: the `/kk` variant must do its per-atom work in a Kokkos
device kernel so it runs on the GPU under a GPU backend. A `/kk` variant that just `sync`s
to host, calls the CPU base-class hook, and marks `modified` is a **fake port** and is
forbidden. The point of porting a thermostat/barostat or forcing fix is **keeping per-atom
data resident on the device** -- a fake port *forces a full `x`/`v`/`f` transfer every
step* and is slower than no `/kk` at all (without one, the integrator can sometimes keep
data on-device across the fix), and it misleads users into thinking the step is
GPU-resident. Nine such fakes (`fix deform/pressure`, `heat`, `indent`, `move`,
`press/berendsen`, `press/langevin`, `restrain`, `temp/csld`, `temp/csvr`) were once
committed and had to be deleted again (the "remove fake KOKKOS ports" commit). If a fix's
per-atom work genuinely cannot run on-device (host-only variable eval, plain `double **`
arrays like `fix move`'s `xoriginal`, external library coupling like `colvars`), **leave
it unported** rather than shipping a delegating wrapper.

## Internal helper computes/fixes must be KOKKOS too

A KOKKOS style that creates an internal helper (a thermostat's temperature, a barostat's
pressure, ...) must ensure the helper is also a KOKKOS style, or it runs on the host and
forces a per-step host/device sync (correct but slow). The trap: `Modify::add_compute` /
`add_fix` only apply the `/kk` suffix under `-sf kk` (`suffix_enable`), so a helper created
with a bare name by a *base-class* constructor is non-kk when the style is requested with
an explicit `/kk` suffix without `-sf kk` (or via `fix_modify temp <non-kk>`). Mitigation,
applied to every such style:

- (A) In the KOKKOS constructor, if the helper exists but is not `kokkosable`, delete and
  recreate it with the `/kk` name, matching the group the base used (`group->names[igroup]`
  or `all`). Guarded, so a no-op when `-sf kk` already promoted it; only for helpers that
  *have* a `/kk` variant -- `pressure` has none. A KOKKOS *compute* recreates its helper in
  `post_constructor()` instead (cf. `compute temp/deform/kk`).
- (C) In `init()`, call `KokkosLMP::warn_nonkokkos_compute(lmp, style, temperature,
  "temperature")` (rank-0 warning, detects via the `kokkosable` member, no-ops if kk or
  null) to catch the residual `fix_modify` case.

Keep every `_kk` helper call `kokkosable`-guarded with a non-kk sync fallback
(`if (helper->kokkosable) helper->foo_kk(); else { sync; foo(); ... }`) so a non-kk helper
is never *wrong*, only slow (the empty base `*_kk()` stubs do nothing). The nh-family is
robust by construction: it hardcodes `temp/kk` via `FixNHKokkos` instead of inheriting the
CPU base's bare helper creation -- prefer that pattern for new styles.

## Documentation checklist (pair, fix, and compute styles alike)

1. New files `src/KOKKOS/<prefix>_<name>_kokkos.{h,cpp}` plus `Install.sh` action lines
   (rule 5). No `src/.gitignore` edit is needed: the `/*_kokkos.{cpp,h}` wildcard rules
   already cover KOKKOS files (only non-`_kokkos` package files need entries).
2. Base class: `virtual allocate()` (rule 2) and destructor `copymode` guard (rule 1).
3. `doc/src/<style>.rst`:
   - Add `.. index:: pair_style <name>/kk` (or `fix <name>/kk`, ...) near the top; index
     entries are grouped by base style, then alphabetical by accelerator suffix (`/kk`
     after `/gpu`, before `/omp`).
   - Add `*<name>/kk*` to the `Accelerator Variants:` line, ordered alphabetically.
   - The `.. include:: accel_styles.rst` directive is REQUIRED in every accelerated-style
     .rst (after a `----------` rule, before Restrictions); `make style_check` flags it.
4. `doc/src/Commands_<type>.rst`: add letter `k` (g=GPU, i=INTEL, k=KOKKOS, o=OPENMP,
   t=OPT) to the style's accelerator string, alphabetically and without commas: `ko`,
   not `ok` or `o,k`.

## Common mistakes (beyond the key rules above)

- Wrong `CoulLongTable` template parameter: use `<0>` with `COUL_FLAG=0`, `<1>` with
  `COUL_FLAG=1`.
- Missing friend declarations: the `pair_compute*` helpers are templated free functions;
  the class must declare them (and the `PairComputeFunctor` instantiations) as friends.
- Using `double` instead of `KK_FLOAT` for floating-point fields in `params_*` structs.
- Missing `// NOLINTNEXTLINE` before a `KOKKOS_INLINE_FUNCTION` (required to suppress
  linter warnings in LAMMPS).
- Missing explicit template instantiation at the bottom of the `.cpp` file
  (`template class <Name>Kokkos<LMPDeviceType>;` plus the `LMPHostType` line under
  `#ifdef LMP_KOKKOS_GPU`).
- Not rebuilding with `make purge` after switching between Make and CMake builds.

## Deep dives

- `.github/dev-docs/kokkos-porting-guide.md` -- pair/fix template skeletons and the full
  lessons-learned collection (custom kernels, wall styles, thermostat/barostat and
  temperature-compute port patterns, device RNG, GPU-only debugging methods).
- `.github/dev-docs/kokkos-porting-backlog.md` -- what remains to port and where open work
  lives; status tables go stale, so trust the tree (`ls src/KOKKOS/`, `git log`) over any
  table.
- Build and test commands: `.github/copilot-instructions.md`; the KOKKOS-specific test
  build block is in the porting guide.
