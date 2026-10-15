# Accelerator Porting Backlog (GPU, KOKKOS, OPENMP)

**STATUS TABLES GO STALE.** Trust `ls src/GPU/ src/KOKKOS/ src/OPENMP/` and `git log`
over this file, and regenerate the pair-style coverage table with
`.github/dev-docs/accelerator-coverage-scan.py` (`--summary`, `--missing omp|gpu|kk`,
`--pkg NAME`, `--yaml`) instead of editing counts by hand.
Last verified: 2026-10-15 (develop at 2b85e289b1).

This file replaces the KOKKOS-only `kokkos-porting-backlog.md` (last refreshed
2026-08-28 after the Tier 1 pass, PR #4989).  Its still-open KOKKOS content is carried
over, condensed, in the KOKKOS section; the OPENMP and GPU sections are new.  The
how-to for each package lives elsewhere: `kokkos-porting-guide.md`, `openmp-porting.md`,
and (to be written with the first GPU port) `gpu-porting-guide.md`.  Scheduled batches
and their progress live outside the tree, in the maintainer's plan files.

## Coverage snapshot

Pair styles (the focus of this backlog), from the scan script on 2026-10-15.  "Base"
includes 7 `DEPRECATED` stubs, `hybrid*` and `zero*`.

| | count |
|---|---|
| base pair styles | 330 |
| with `/omp` | 204 |
| with `/kk` | 138 |
| with `/gpu` | 72 |
| with `/intel` / `/opt` | 19 / 16 |
| no accelerator variant at all | 100 (CG-DNA 23, SPIN 7, APIP 7, MACHDYN 4, ML 6, SPH 3, PERI 2, ...) |
| with a `unittest/force-styles` YAML reference | 246 |

Other style kinds under KOKKOS (counted 2026-08-28; the `/kk` column counts source
files): Fix 343/87, Compute 188/25, KSpace 27/2, AtomVec 30/12, bonded 97/71, Region
10/7, Minimize 10/4.  `NPair`/`NStencil`/`NBin` are not a gap (see "Neighbor" below).

## Exclusions

**Pending pull requests (do not touch; check `gh pr list` before starting anything):**
`pace` (/omp, #5243), `granular` (/kk, #5208 and #5187), the oxDNA/oxRNA family (/kk,
#5141), `dispersion/d3` (/omp and /kk, #5056), MTP (#5023), `reaxff` and `qeq/reaxff`
(#5242), `lj/class2/*soft/gapsys` (#5162), `lubricate/poly` and `brownian/poly`
(#5040), `lj/cut/coul/long2` and the `pair_kokkos.h` neighbor-cluster changes (#5086).

**Standing decisions:** PERI is being replaced by BPM (no new work on `peri/*`); SPIN
(`fix nve/spin` is a sequential sweep, a threaded or device `compute()` gains <= 1.2x);
the lubricate family is slated for removal; AIREBO/REBO/rebomos/COMB/COMB3/BOP/LCBOP are
not ported by any accelerator package and will not be; ML-QUIP is keep-building-only;
external-library wrappers (`lepton`, `python`, `kim`, `quip`, `hdnnp`, `mbx`, `mgpt`,
`runner`) have nothing to offload; MSM pair styles are pointless for KOKKOS and GPU
without `msm/kk` or `msm/gpu`; no new `/intel` variants (the core developers only
maintain INTEL; contributor INTEL code is welcome).

## Effort scale

Used by every table below.  Estimates include the docs and the unit-test run.

| level | meaning |
|---|---|
| TRIVIAL | under 2 hours: near-mechanical clone of a sibling variant |
| SMALL | half a day: clone plus a changed formula, a second RNG, a torque term, ... |
| MEDIUM | 1-2 days: new kernel shape, per-thread or device scratch, two-phase comm |
| LARGE | several days, or blocked on a framework piece listed next to it |

Per-package cost floor, before any physics: OPENMP = 2 files (`src/OPENMP/
pair_X_omp.{h,cpp}`, registration automatic); KOKKOS = 2 files + `src/KOKKOS/Install.sh`;
GPU = 6 files (`src/GPU/pair_X_gpu.{h,cpp}`, `lib/gpu/lal_X.{h,cpp,cu}`,
`lib/gpu/lal_X_ext.cpp`) + a declaration block in `lib/gpu/include/lammps_gpu.h` +
`src/GPU/Install.sh` + `src/Purge.list`.  All three: the `(o)`/`(k)`/`(g)` letter in
`doc/src/Commands_pair.rst`, an `.. index::` line and the "Accelerator Variants" entry
of the pair page, no `versionadded`.

Validation is free when the base style has a YAML reference: `test_pair_style` runs
the `omp`, `kokkos_omp`, `kokkos_serial`, `kokkos_gpu`, `gpu`, `intel` and `opt` cases
off the CPU reference automatically.  A missing variant SKIPS SILENTLY, so confirm in
`ctest -V` output that the new case ran, and drop stale accelerator tokens from
`skip_tests:` when a port makes them runnable.

## OPENMP

### Framework item: per-thread reduction beyond f and torque

`FixOMP::pre_force` hands every thread a zeroed slice of `atom->f`, `torque`,
`erforce`, `desph` and `drho` (`fix_omp.cpp`, `ThrData::init_force`), and `Atom`
allocates `nthreads` copies of those plus `edpd_flux`, `cc_flux`, `fm`, `fm_long`
(`atom.cpp`, `threadflag`).  But `ThrOMP::reduce_thr` (`thr_omp.cpp`) folds back only
`f` and `torque`, so no `/omp` pair style has ever used the other slices.  The work:
add `data_reduce_thr` calls for `erforce`, `desph`, `drho`, `edpd_flux` and `cc_flux`
(stride `cc_species`) in the same `last_omp_style` block as `torque`, extend
`ThrData::init_force` / `FixOMP::pre_force` with `edpd_flux` and `cc_flux` slices, add
`get_erforce()`/`get_edpd_flux()`/`get_cc_flux()`.  Not in the `last_pair_hybrid`
early-reduction branch (that exists only for `PairHybrid`'s `virial_fdotr_compute`).
Consumers accumulate into `thr->get_drho()` etc. exactly as the base does and never
zero or reduce the slice themselves (zeroing "its" slice at `eval()` start, the
`init_eam` idiom, erases the previous hybrid sub-style's contribution).  Ordering is
already right: the fold completes before `comm->reverse_comm()`, the SPH/EDPD/EFF
integrators and `Min::energy_force -> min_xf_get`.  Known gap shared with torque:
`respa/omp` reduces only `f`.  The YAML harness compares only `f`: validate with
`-pk omp 1` vs `-pk omp 4` vs the plain style on `examples/PACKAGES/sph`, `dpd-meso`
and `eff` decks, and a negative control (remove one reduction, watch it diverge).
EAM/CD scratch (`rhoB`, `D_values`) stays style-owned like `pair_eam_omp`; a
`ThrData::init_cdeam()` next to `init_eim()` is the natural hook.

### Pair styles without /omp

| effort | styles (template to clone) |
|---|---|
| TRIVIAL | `lj/cut/soft/gapsys` (coul/cut/soft/gapsys/omp + lj/cut/soft/omp); `nb3b/screened` (derive from `PairNb3bHarmonicOMP`, set the variant, override `threebody()`: `eval()` calls the virtual); `sw/angle/table` (copy the sw/omp eval, call `threebody_table`; diamond inheritance rules out deriving from `PairSWOMP`); `lj/long/coul/long/dielectric` (full list, writes only atom i; lj/cut/coul/long/dielectric/omp + the dispersion branch of lj/long/coul/long/omp); `bpm/peri` (bpm/spring/omp; `vfrac` is read-only); `sph/heatconduction` (desph only); `sph/rhosum`, `mdpd/rhosum` (full list, `rho[i]` only, then forward comm); `smatb/single` (after smatb) |
| SMALL | `dpd/coul/slater/long` (dpd/omp per-thread RNG + coul/slater/long/omp); `lj/long/dipole/long` (lj/long/coul/long/omp + lj/cut/dipole/long/omp: torque, `ev_tally_xyz`); `edip/multi` (edip/omp; scratch already sized `nthreads*5*64`); `meam/sw/spline` (meam/spline/omp two-phase; 19 example decks); `smatb` (eam/omp two-phase with reverse and forward comm); `multi/lucy` (eam/omp phases + table/omp lookup); `sph/idealgas`, `sph/taitwater`, `sph/taitwater/morris`, `sph/lj` (lj/cut/omp + `get_drho()`/`get_de()`); `sdpd/taitwater/isothermal` (dpd/omp RNG + drho); `edpd`, `tdpd` (dpd/omp + flux slices); `multi/lucy/rx`, `table/rx` (pair-owned per-thread `uCG`/`uCGnew`, reduced before `fix eos/table/rx` reverse-comms them); `uf3` (sw/omp; `neighshort` per thread) |
| MEDIUM | `eam/cd` + `eam/cd/old` (CDEAM scratch, two or four comm stages inside `eval()`); `polymorphic` (tersoff/omp; member scratch per thread); `kolmogorov/crespi/full` (drip/omp: serial ILP neighbor prep, threaded eval); `exp6/rx` (1.3k lines, per-atom parameter arrays); `eff/cut` (new `ev_tally_eff_thr`, `erforce` fold before `min_xf_get()`, and it has no YAML reference: needs an `atom_style electron` deck first); `snap` (one `SNA` instance per thread; `/kk` shows the decomposition) |
| LARGE / not worth it | `bop` (serial `gneigh` + 1200 lines of `SigmaBo`/`PiBo` on member scratch writing `atom->f`), `comb3`, `mliap` (Python models, GIL), `pod`, `amoeba`/`hippo` (2.5k lines, comm + map + fix), `mesocnt`(+`viscous`), `smd/*`, `surf/granular`, `granular/superellipsoid`, `rheo`, `lubricateU`(+`poly`), APIP styles, `meam` (library), `spin/*` (SMALL each, <= 1.2x gain), `granular` (wait for the #5208 kernel refactor) |

Orphans and drift: none found 2026-10-15 (every `pair_*_omp.cpp` has a base; the 30
`/omp` files untouched since the 2022 bulk reformat whose bases changed in 2025-2026
were spot-checked, all cosmetic).  Maintenance method and bug classes:
`openmp-porting.md`.

## KOKKOS

### Framework facts (carried over, still true)

- A device-callable atom map exists: `AtomKokkos::map_kokkos<DeviceType>()`
  (`fix_shake_kokkos`, `fix_cmap_kokkos`, `fix_rigid_small_kokkos` use it).
- Region + atom-style variables are solved in `fix_setforce_kokkos.cpp`
  (`dynamic_cast<KokkosBase *>` + `match_all_kokkos`, host-evaluated variables in a
  DualView); the region type list lives once in `region_kokkos_styles.h`.
- Group reductions exist on device (`group_kokkos.h`: mass, xcm, vcm, angmom, inertia,
  gyration, gyration tensor).
- Neighbor stencils are deliberately host-computed and memcpy'd.
- `PPPMKokkos` is ik-differentiation only (`kspace_modify diff ad` errors).
- Per-atom density with device forward+reverse comm exists in `pair_eam_kokkos`
  (`pack/unpack_forward_comm_kokkos`, `pack/unpack_reverse_comm_kokkos`): the
  template for every EAM-shaped style below.
- Missing: a uniform device per-atom virial/energy accessor across the force
  hierarchies; device reverse-comm for computes; KOKKOS `fix STORE/ATOM`; a KOKKOS
  local-data (`array_local`) framework; `compute chunk/atom` on device.

### Pair styles without /kk

| effort | styles (device template to clone) |
|---|---|
| SMALL (`pair_kokkos.h` template) | FEP leftovers `coul/cut/soft/gapsys`, `lj/cut/soft/gapsys`, `lj/class2/coul/cut/soft`, `lj/class2/coul/long/soft` (coul/cut/soft/kk, lj/class2/soft/kk, lj/cut/coul/long/soft/kk); `coul/esp`, `lj/cut/coul/esp` (coul/long/kk, lj/cut/coul/long/kk including the table branch); `dpd/coul/slater/long` (dpd/kk device RNG + coul/slater/long/kk; YAML needs the `kokkos_*_devicerng` tokens); `lj/sf/dipole/sf`, `mesomem/dipole` (lj/cut/dipole/cut/kk hand kernel with torque); `tip4p/long/soft`, `lj/cut/tip4p/long/soft` (lj/cut/tip4p/long/kk); `coul/ctip` (pairwise + self terms; the `qeq_*` parameter arrays only feed `fix qeq/ctip` on the host) |
| MEDIUM (many-body, template exists) | `nb3b/harmonic` + `nb3b/screened` (copy sw/kk, delete twobody, switch on the variant); `gw` + `gw/zbl` (tersoff/kk); `vashishta/table` (vashishta/kk + device table); `sw/angle/table`, `threebody/table` (sw/kk + device table views as in table/kk); `edip` (sw-shaped full list); `tersoff/table`; `extep`; `atm` (clean triple loop over one list; device `ev_tally3` exists in sw/vashishta/uf3 kk); `kolmogorov/crespi/z`, `lebedeva/z` (hand kernels, anisotropic force, no normals); `coul/exclude` (hand kernel: no `cutsq` guard, skips `sbmask(j)==0`); `coul/streitz` (full list; `fix qeq/slater` stays host) |
| LARGE (framework piece first) | EAM-shaped density styles `eam/cd`, `local/density`, `eim`, `meam/spline`, `meam/sw/spline`, `smatb`(+`single`), `multi/lucy` (MEDIUM each once the first reuses the eam/kk comm pattern); DRUDE `thole`, `lj/cut/thole/long`, `coul/tt` (partner indices precomputed per reneighbor through the device map); `hbond/dreiding/*` (topology walk via `tagprev`); ASPHERE `gayberne` then `resquared` (`AtomVecEllipsoidKokkos` + `MathExtraKokkos::quat_to_mat` exist; ylz/kk is the only ellipsoid pair so far); INTERLAYER `ilp/graphene/hbn`, `ilp/tmd`, `saip/metal`, `saip/metal/tmd`, `aip/water/2dm`, `kolmogorov/crespi/full`, `drip` (all gated on `ILP_neigh()` + `calc_normal()` on device; once the first exists the rest ride along); `lj/long/coul/long`, `buck/long/coul/long`, `lj/long/tip4p/long`, `lj/long/dipole/long`, `lj/cut/dipole/long` (pointless before `pppm/disp/kk` or a dipole kspace on device); DIELECTRIC (7 styles: need `atom_vec_dielectric_kokkos` and `fix polarize/*`); atom-vec gated: SPH (6 + `sdpd/taitwater/isothermal`), DPD-MESO (`edpd`, `mdpd`, `mdpd/rhosum`, `tdpd`), EFF, MACHDYN, RHEO (`rheo`, `rheo/solid`), BODY, `line/lj`, `tri/lj`, MESONT, BPM (`bpm/spring`, `bpm/peri`: `atom_vec_bpm_sphere_kokkos`); `polymorphic`, `smtbq`, `e3b`, `brownian/poly` (after #5040) |
| do not port | MSM pair styles (8), SPIN (7), AIREBO family, `rebo`, `rebomos`, `comb`, `comb3`, `bop`, `lcbop`, lubricate family, external wrappers, `srp`, `list`, `tracker`, `dsmc`, `zero` |

Atom styles (`sph`, `peri`, `edpd`/`mdpd`/`tdpd`, `electron`, `smd`, `apip`, `oxdna`)
remain the highest-leverage *gating* items: `AtomKokkos::new_avec` hard-errors on them,
and each needs new masks in `atom_masks.h`.  Do **not** port `atom_style template`
(`neighbor_kokkos.cpp` errors on molecule templates because the KK neighbor build reads
special bonds from the per-atom `special` list only; porting the vec alone would make
neighbor lists silently wrong).  `bond special` is not portable (needs a device
`Pair::single()`).

### Tier 0 -- framework work that unblocks the most

1. **`NBinMulti` on device** + a `MULTI` template parameter on `NPairKokkos`.
   `neighbor multi` is a hard error under KOKKOS, which costs polydisperse granular and
   colloidal users the O(N) optimization exactly where `gran/*/kk`, `colloid/kk`,
   `yukawa/colloid/kk` and `wall/gran/kk` already run.  Stencils come free (host
   `NStencilMulti` via the existing memcpy).  Highest-value neighbor work by far.
2. **KOKKOS `fix STORE/ATOM`** (DualView `astore` + `pack_exchange_kokkos`): unblocks
   `compute displace/atom` and `compute msd`.
3. **Device reverse-comm for computes** (`KokkosBase` + `CommKokkos`): unblocks
   `compute contact/atom`, `snad/atom`, `snav/atom`.
4. **Uniform device per-atom virial/energy accessor**: unblocks `compute stress/atom`,
   `pe/atom`, `centroid/stress/atom`, `heat/flux` and the meta-computes (`reduce`,
   `slice`, `global/atom`, `fix pair`).
5. **Lift the `kspace_modify diff ad` restriction** in `pppm_kokkos.cpp`.
6. **Region dispatch refactor** -- partly done (type list in `region_kokkos_styles.h`);
   replacing the `dynamic_cast` chain with a virtual device interface is what is left.

### Non-pair styles: Tier 2 (moderate effort, clear payoff)

- **QEQ family (6 fixes + `qeq/comb`)**: CSR `H` over the neighbor list + CG with a
  forward comm per iteration, exactly what `fix_qeq_reaxff_kokkos.cpp` does; the shared
  `FixQEq` base means one port covers all.  The largest block of unblocked device work.
- `fix indent` (four per-atom loops; ~1140 lines now, no longer small);
  `nvt/npt/nph/asphere` as `FixNHAsphereKokkos : FixNHKokkos` + `compute
  temp/asphere/kk`; `pppm/stagger/kk` (run the existing kernels twice with a shifted
  `shiftone`: the most realistic KSPACE target); `fix msst`, `fix rattle` (subclass of
  the ported `FixShake`), `fix wall/table`, `wall/ees`, `wall/reflect/stochastic`,
  `fix gld`, `fix ffl`; `compute rdf` (histogram with `atomic_add`), `compute sna/atom`
  (`sna_kokkos_impl.h` has the descriptor on device).

### Non-pair styles: Tier 3 (hard, moves flagship workloads)

- **`pppm/disp`** (8860 lines): realistic scoping is `function[0]` (coulomb) +
  `function[1]` (geometric mixing) only, hard-error on arithmetic/none; even that is
  ~2x `pppm_kokkos.cpp`.  Prerequisite for the `lj/long`/`buck/long` pair styles.
- `fix ttm`, `ttm/grid`; `ewald/kk` (O(N^1.5) but the easiest kspace: two
  `parallel_reduce` over kmax x nlocal plus a force kernel; has `/gpu` and `/omp`).

### Do NOT port (KOKKOS)

- No device work to offer: `compute pe`, `compute pressure` (no per-atom loop), `fix
  ave/time`, `fix vector`, `fix print`, `fix set`, `fix store/state`, `compute
  property/atom`, all `Command` and `Dump` styles.
- Meta-styles over other styles' host pointers (a `/kk` would be a fake port):
  `compute reduce`, `reduce/region`, `slice`, `global/atom`, `fix pair`, `fix ave/atom`.
  They fall out of Tier 0 item 4.
- Atom count or topology changes: `fix deposit`, `evaporate`, `pour`, `append/atoms`,
  all of MC; `fix tfmc`.
- Inherently serial: `fix thermal/conductivity`, `fix viscosity`, `compute
  born/matrix`, `stress/mop`, `fabric` (host `pair->single()`), `fix numdiff`,
  `numdiff/virial`, `fix gle`.  `fix move`, `fix deform/pressure` (host scalar work);
  `fix restrain` (portable, no payoff).
- Whole TALLY package (`Pair::add_tally_callback` host callback per pair).
- The whole SPIN package (sequential symplectic sweep with `forward_comm` inside the
  loop; all six `pair_spin_*` unported; `min spin*` likewise).
- `kspace msm`, `msm/cg` and every `*/msm` pair style; ELECTRODE (dense capacitance
  matrix on the host); `atom_style body`; `fix rigid` (big) and `rigid/nve|nvt|npt|nph`
  (`rigid/small` is ported and recommended); `pair srp`; `fix ilves/omp`-style
  per-cluster threading (benchmarked and rejected); r-RESPA (four independent
  blockers, recommend against without a concrete user).
- `fix npt/kk` and `nph/kk` are complete full ports; nothing remains on Nose-Hoover.

### Neighbor

Do not track `NStencilStyle` or `NPairStyle` as a gap: stencils are host-computed by
design, ~half the non-kk `NPairStyle` registrations are `/omp` variants, the
`atomonly` micro-variants fold into the templated `NPairKokkos`, `nsq` is deliberately
rejected, `multi/old` is gone, and respa lists cannot be requested.  The honest gap is
`multi` binning (Tier 0 item 1) plus `skip/size/off2on` and `skip/size/off2on/oneside`
(granular hybrid setups, `fix wall/gran/region`, GRANSURF) as a distant second.

### History

Tier 1 (PR #4989, 2026-08/09): 23 pair styles on the `pair_kokkos.h` template (8 FEP
soft-core, 8 CORESHELL `/cs`, 7 singletons), 4 many-body derivations (`sw/mod`,
`tersoff/mod/c`, `gran/hooke`, `gran/hertz/history`), 4 bonded table styles, 13 fixes,
7 computes, 5 regions, 2 minimizers.  Lessons worth keeping: the `/cs` B0..B5 erfc
series pairs with its own EWALD_P (9.95473818e-1), not `EwaldConst::EWALD_P`; the
`EPS_EWALD_SQR` correction reaches the LJ term on the CPU through a shared `r2inv` but
only the Coulomb term in the template (eight YAML epsilons were raised ~4x the measured
7.7e-12); device RNG streams cannot reproduce the CPU Marsaglia stream (YAML skips
`kokkos_*` for `baoab`, `brownian`, `langevin`, `gjf`).  CPU bugs fixed by that pass:
`lj/smooth/linear` (+`/omp`, `lj/sf`) tallied the pair energy without `factor_lj`;
`~PairLJCharmmfswCoulCharmmfsh` restored `qqr2e` before checking `copymode`.  Left
alone and reproduced by `/kk`: `morse/soft` does not subtract `offset` in `compute()`
while `single()` does; `compute entropy/atom`'s first bin has `rbinsq[0] = 0`.
Verification decks for all Tier 1 styles: `kokkos-tier1-verification/`.

## GPU

### Template facts

- One style = `src/GPU/pair_X_gpu.{h,cpp}` (subclass of the CPU class overriding ctor,
  dtor, `compute`, `init_style`, `memory_usage`, plus `reinit` when the base has it),
  `lib/gpu/lal_X.{h,cpp,cu}` (device class + kernels `k_X` / `k_X_fast`),
  `lib/gpu/lal_X_ext.cpp` (C entry points `X_gpu_init/compute/compute_n/clear/bytes` in
  `namespace LAMMPS_GPU`), a declaration block in `lib/gpu/include/lammps_gpu.h`
  (included on both sides, so signatures are checked), two `action` lines in
  `src/GPU/Install.sh`, an entry in `src/Purge.list`.  CMake globs `lib/gpu/*.cpp` and
  `*.cu` for CUDA, OpenCL and HIP (`cmake/Modules/Packages/GPU.cmake`); the kernel
  symbol is the file name minus `lal_`; the host `.cpp` includes `X_cl.h` /
  `X_cubin.h`.  The legacy `lib/gpu/Makefile.*` are gone.  Clean template commit:
  086cf49a8c (`coul/slater/long`, 2023-12), minus three API changes since
  (`lammps_gpu.h` consolidation 6de3baf5bf/e8059a06dd, `host_start`/`cpu_time` removal
  21d55a1540, `periodicity` argument baf16f9117): copy the structure from the commit
  and the code from the current `lal_coul_slater_long.*`.
- Base families in `lib/gpu/lal_base_*`: `atomic` (lj, born, buck, morse, mie, eam,
  table, zbl, yukawa), `charge` (q array + texture, `threads_per_charge`; `qqrd2e` and
  `g_ewald` are per style), `dipole`, `dpd` (uploads `v` + `tag`, in-kernel Saru RNG
  seeded from the tag pair), `sph`, `three` (2-body + three-center + three-end kernels,
  `REQ_FULL|REQ_GHOST` lists), `ellipsoid`, `amoeba`.
- Long-range Coulomb is the analytic erfc polynomial only; no kernel reads
  `ncoultablebits`, so every `*_table` YAML lists `gpu` in `skip_tests` and
  `init_tables` is called on the host only for the inherited `single()`.
- Both newton settings work (full lists, forces on local atoms only; the per-style
  newton checks went in 2021).  `eatom`/`vatom` are accumulated on device and added in
  `FixGPU::post_force`.  Hybrid works with host neighbor lists (`neigh no`) only.
  `respa_enable = 0` everywhere.  Triclinic works; `neigh_modify exclude` and
  `atom_style template` do not with device neighbor builds.
- Precision at configure time: `-D GPU_PREC=single|mixed|double` (default mixed).
  The YAML harness relaxes `epsilon` x7.5 (double), x5e8 (mixed), x1e10 (single):
  validate new kernels against a double-precision build, then run mixed and HIP.
- `src/GPU/Install.sh` is missing the `action` lines for nine existing styles
  (`coul/slater/long`, the SPH styles, the soft styles among them); harmless for CMake
  builds, fix when touching the file.

### Pair styles without /gpu

| effort | styles (closest `lal_` kernel) |
|---|---|
| SMALL (charge family) | `lj/charmmfsw/coul/long` (lal_charmm_long + force switching), `lj/charmmfsw/coul/charmmfsh`, `lj/charmm/coul/charmm/implicit` (lal_charmm); `lj/class2/coul/cut` (lal_lj_class2_long + lal_lj_coul); `coul/wolf`, `lj/cut/coul/wolf` (lal_born_coul_wolf); `born/coul/dsf`, `born/coul/dsf/cs` (lal_coul_dsf + lal_born); `coul/slater/cut` (lal_coul_slater_long); FEP `coul/cut/soft`, `coul/long/soft`, `lj/cut/soft`, `lj/class2/soft`, `morse/soft`, `lj/charmm/coul/long/soft` (lal_lj_coul_soft, lal_lj_coul_long_soft); CORESHELL `buck/coul/long/cs`, `lj/cut/coul/long/cs`, `lj/class2/coul/long/cs`, `coul/wolf/cs` (lal_born_coul_long_cs, lal_coul_long_cs, lal_born_coul_wolf_cs); `lj/gromacs/coul/gromacs` (lal_lj_gromacs); `nm/cut/coul/cut`, `nm/cut/coul/long` (lal_lj_coul, lal_lj_coul_long); `coul/diel`, `coul/shield`, `coul/cut/global` (lal_coul); `lj/switch3/coulgauss/long`, `mm3/switch3/coulgauss/long`, `coul/esp`, `lj/cut/coul/esp` (lal_lj_coul_long); `pedone`, `buck6d/coul/gauss/dsf`, `buck6d/coul/gauss/long`, `lj/class2/coul/cut/soft`, `lj/class2/coul/long/soft`, `coul/cut/soft/gapsys` |
| SMALL (atomic family) | `nm/cut`, `nm/cut/split`, `harmonic/cut`, `wf/cut`, `momb`, `lj/pirani`, `lj/mdf`, `buck/mdf`, `lennard/mdf`, `lj/smooth/linear` (+`lj/sf` alias), `morse/smooth/linear`, `lj/relres`, `cosine/squared`, `born/gauss`, `gauss/cut`, `lj/cut/soft/gapsys` (lal_lj and siblings); with a radius upload (lal_colloid / lal_yukawa_colloid pattern): `lj/cut/sphere`, `lj/expand/sphere`; DPD `dpd/ext`, `dpd/ext/tstat` (lal_dpd) |
| MEDIUM | `tip4p/long`, `tip4p/cut`, `lj/cut/tip4p/cut` (lal_lj_tip4p_long M-site kernel); `mesomem/dipole` (lal_dipole_lj); three family `sw/mod` (lal_sw), `tersoff/mod/c` (lal_tersoff_mod), `vashishta/table` (lal_vashishta), `nb3b/harmonic`(+`screened`), `gw`(+`zbl`) (lal_tersoff), `sw/angle/table`, `threebody/table`, `edip`; `ylz` (lal_gayberne); `lj/long/coul/long` (lal_lj_coul_long + the analytic dispersion terms; the kspace stays on the host anyway); density styles `eam/cd`, `local/density`, `eim` (lal_eam) |
| LARGE / not now | `atm` (no triple-loop template), `meam`, ML (`snap`, `pace`, `pod`, `uf3`, `mliap`), the ILP family, `granular` (owner: the GPU package maintainer), BPM, SPIN, oxDNA, DIELECTRIC, `eff/cut` (requested in #3325; electron radius data on device), `amoeba`/`hippo` already exist |

Known gaps in existing `/gpu` variants (from `skip_tests`): no Coulomb tables (by
design), `tersoff*` reject `shift`, `sw` rejects `threebody off`, hybrid and
hybrid/scaled need host neighbor lists, the DPD and QEQ families skip because of RNG
or host-side solvers; `edpd/gpu` and `mdpd/gpu` are unverified since the SPH GPU
crash fixes of #5212.

## Tests

- Mechanics (`unittest/force-styles/test_pair_style.cpp`): the skip token in
  `skip_tests:` is the gtest case name (`omp`, `gpu`, `intel`, `opt`, `kokkos_omp`,
  `kokkos_serial`, `kokkos_gpu`, plus `<case>_double|_mixed|_single`,
  `<case>_devicerng`, `<case>_full`, `single_thread`, `host_neigh`, `gpu_no_single`,
  `gpu_no_mixed`, `single`, `extract`, `nofdotr`, `vatom_only`, `numdiff`).  KOKKOS
  cases run `-k on t 4 -sf kk` (OpenMP) / `-k on t 1` (Serial) / `-k on g 1 -pk kokkos
  neigh half newton on` (device), epsilon x5, then x2e9 mixed / x1e10 single.  GPU runs
  `-sf gpu` (`-pk gpu 0 neigh no` for hybrids), newton off only.  OPENMP runs
  `-pk omp 4 -sf omp`, epsilon x5.  Generation: `test_pair_style <template.yaml> -g
  <new.yaml>` (first argument is a hand-edited copy of a related YAML), remove the
  `generated` tag after review, re-run cmake.
- Pair styles WITHOUT a YAML reference and what a reference needs: `eff/cut` (new
  deck + data file: `atom_style electron`, `units electron`; the harness compares `f`
  only, not `erforce`); `line/lj`, `tri/lj` (new 2d line / tri decks, or `in.empty` +
  `set ... line|tri`); `amoeba`, `hippo` (new deck: `atom_style amoeba`, `fix
  property/atom`, `.prm`/`.key` files; newton on, 3d); `body/nparticle` (new data file
  with a Bodies section; existing body decks are rounded/polygon and polyhedron);
  `mesocnt` (new deck: CNT bead chains, `special_bonds lj 1 1 1`, newton on + full
  list); `rheo` (`in.rheo` exists; needs the `fix rheo`, `rheo/pressure`,
  `rheo/viscosity` stack in `post_commands`; `fix rheo` integrates on its own);
  `surf/granular` (line/tri with radius, newton on, ghost velocities);
  `granular/superellipsoid` (`in.empty` + `atom_style ellipsoid superellipsoid`;
  already covered by `unittest/granular/`); `smd/*` (own integrator fix + `tlsph`
  group: conflicts with the harness, hardest); `peri/*` (covered by
  `test_pair_peri.cpp` on purpose); `mesont/tpm` (deprecated stub).  Any YAML built on
  `in.empty` must skip `single`.
- The YAML harness compares only `f`, energy and stress: per-atom accumulators
  (`drho`, `desph`, `erforce`, `edpd_flux`, `cc_flux`, `efield`) need an OutputStyle
  YAML (`test_output_style.cpp`, e.g. clones of `compute-sph_rho_atom.yaml` with the
  pair style in `post_commands`) or an example-deck comparison.

## Bugs found by the 2026-10-15 survey

1. `src/OPENMP/pair_edip_omp.cpp`: the per-thread `preForceCoord` slice was offset by
   `tid * 64` while each thread uses `5 * 64` entries (the base allocates
   `5 * nthreads * 64`), so threads overlapped once an atom had more than 12 neighbors
   in the cutoff-function range.  `edip`, `edip/omp` and `edip/multi` also never
   checked the neighbor count against the fixed 64-entry scratch.
2. `src/DPD-MESO/pair_edpd.cpp`: `read_restart_settings()` seeded `randomT` with the
   same seed as `random` (`settings()` uses `2*seed`), so after a restart the two noise
   streams were identical.
3. `src/MANYBODY/pair_nb3b_harmonic.cpp`: substyle string typo `"nb3n/harmonic"`.
4. Design question, not a bug fix: `pair bpm/peri` tallies `fpair*vfrac_eff` into the
   virial on purpose but does not set `no_virial_fdotr_compute`, so the global virial
   differs between the fdotr path and the per-pair tally path.

Items 1-3 are fixed on their own branch; see the commit history.
