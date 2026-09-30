# KOKKOS Porting Backlog

**STATUS TABLES GO STALE.** Trust `ls src/KOKKOS/` and `git log` over this file.
Last verified: 2026-08-28 (end of the Tier 1 pass).

Derived from an exhaustive scan: every style-registration macro (20 of them --
`PairStyle`, `FixStyle`, `ComputeStyle`, `KSpaceStyle`, `RegionStyle`, ...) across
all `src/**/*.h`, minus everything carrying a `/kk` registration, scored on three
signals -- an existing `/gpu`/`/intel`/`/omp`/`/opt` variant (someone already found
the style hot), occurrence count across `examples/` + `bench/`, and whether a
`unittest/force-styles/tests/*.yaml` reference exists (96 unported pair styles and
18 unported kspace styles have one, so a port is cheap to validate).

## Coverage snapshot

Counts are as of the end of the Tier 1 pass; the `/kk` column is
`grep -l "^<Kind>Style(" src/KOKKOS/*_kokkos.h | wc -l`, which counts source
files, not registrations.

| Kind | base | `/kk` | real gap\* |
|---|---|---|---|
| Fix | 343 | 87 | ~213 |
| Pair | 328 | 135 | ~181 |
| Compute | 188 | 25 | ~149 |
| KSpace | 27 | 2 | 23 |
| AtomVec | 30 | 12 | 18 |
| Bond / Angle / Dihedral / Improper | 97 | 71 | 18 |
| Region | 10 | 7 | 2 |
| Minimize | 10 | 4 | 4 |
| NPair + NStencil + NBin | 103 | 19 | **1** (see "Neighbor" below) |
| Command, Dump | 95 | 2 | 0 -- setup-time / host I/O |

\* excludes `DEPRECATED` aliases, `zero` styles, upper-case internal styles, and
external-library wrappers.

History: pair Groups 1-4 (~30 styles), fix Groups A-C, the TIP4P family (`tip4p/cut`,
`tip4p/long`, `lj/cut/tip4p/cut`, `lj/cut/tip4p/long`, `pppm/tip4p`), fix Groups D-E
(`spring/rg`, `ti/spring`, `addtorque/group`, `gjf`, `press/berendsen`,
`press/langevin`, `temp/csvr`, `temp/csld`), `fix heat`, and computes
`temp/partial`, `temp/profile`, `temp/ramp`, `ke/atom` are all **merged**.  Nine fake
(sync-to-host) fix ports were committed and deleted again; see the "remove fake
KOKKOS ports" commit.

## Framework facts worth knowing before you scope a port

These are easy to get wrong from memory, and several older backlog entries were
written before they were true:

- **A device-callable atom map exists.** `AtomKokkos::map_kokkos<DeviceType>()`
  (`src/KOKKOS/atom_kokkos.h:144`, `atom_map_kokkos.cpp`), already used by
  `fix_shake_kokkos`, `fix_cmap_kokkos` and `fix_rigid_small_kokkos`.  "It calls
  `atom->map()`" is a cost, not a blocker.
- **Region + atom-style variables are solved** in `fix_setforce_kokkos.cpp`: it
  errors out via `dynamic_cast<KokkosBase *>` on an unported region, calls
  `regionKKBase->match_all_kokkos(groupbit, k_match)`, and host-evaluates atom-style
  variables into a DualView the kernel reads.  Any `setforce`-shaped fix inherits
  the recipe for free.
- **Group reductions exist on device.** `src/KOKKOS/group_kokkos.h` provides
  `mass_kk`, `xcm_kk`, `vcm_kk`, `angmom_kk`, `inertia_kk`, `gyration_kk` and
  `gyration_tensor_kk`.
- **Neighbor stencils are deliberately host-computed.**
  `npair_kokkos.cpp::copy_stencil_info()` calls the base `NPair::copy_stencil_info()`
  and memcpy's the result into device views.  A stencil is a handful of integer
  offsets computed once per re-bin; there is nothing to parallelize.
- **`PPPMKokkos` is ik-differentiation only** -- `pppm_kokkos.cpp:145` errors on
  `kspace_modify diff ad`.  This shapes every PPPM-adjacent estimate.
- **Missing framework pieces** that gate whole clusters: no uniform device accessor
  for per-atom virial/energy across the seven force hierarchies; no device
  reverse-comm for computes; no KOKKOS `fix STORE/ATOM`; no KOKKOS local-data
  (`array_local`) framework; no `compute chunk/atom` on device.

## Tier 0 -- framework work that unblocks the most

1. **`NBinMulti` on device** + a `MULTI` template parameter on `NPairKokkos`.
   `neighbor multi` is a hard error under KOKKOS (`neighbor_kokkos.cpp:125`), which
   costs polydisperse granular and colloidal users the O(N) collection optimization
   -- exactly the workloads where `gran/hooke/history/kk`, `colloid/kk`,
   `yukawa/colloid/kk` and `wall/gran/kk` already run on device.  The stencil side
   is free (see above): reuse the host `NStencilMulti` stencils via the existing
   memcpy.  **Highest-value neighbor work by a wide margin.**
2. **KOKKOS `fix STORE/ATOM`** (DualView `astore` + `KokkosBase::pack_exchange_kokkos`,
   `fix_spring_self_kokkos` pattern).  Unblocks `compute displace/atom` and
   `compute msd` -- that is the right unit of work, not the computes.
3. **Device reverse-comm for computes** (`KokkosBase` + `CommKokkos`).  No compute
   anywhere does `pack/unpack_reverse_comm_kokkos` today.  Unblocks
   `compute contact/atom`, `snad/atom`, `snav/atom`.
4. **Uniform device per-atom virial/energy accessor** (`k_vatom`/`k_eatom` live
   privately inside each KOKKOS subclass).  Unblocks `compute stress/atom`,
   `pe/atom`, `centroid/stress/atom`, `heat/flux`, and the meta-computes
   (`reduce`, `slice`, `global/atom`, `fix pair`) that read other styles'
   `vector_atom`/`array_atom` host pointers.
5. **Lift the `kspace_modify diff ad` restriction** in `pppm_kokkos.cpp:145`.  Not a
   new style, but `diff ad` is the memory-efficient recommendation for large PPPM
   runs -- arguably higher value than any new kspace style.
6. **Region dispatch refactor -- partly done.** The concrete-type list that
   `fix wall/region/kk` templates its functor on now lives once, in
   `region_kokkos_styles.h`, and drives both the dispatch and the error message, so
   adding a region/kk style is a one-line change there.  The dispatch is still a
   `dynamic_cast` chain with one template instantiation per region type times
   device/host times precision; replacing it with a virtual device interface is
   what is left.

## Tier 1 -- mechanical, template already in tree

**Complete** on branch `more-kokkos-porting`, except for the two items listed
under "Not done" at the end.  Everything below was verified either against the
existing `unittest/force-styles` YAML references (which exercise every `/kk`
variant automatically) or, where no reference exists, against the CPU style on a
hand-written deck run serial and under `mpirun -np 2`.

### Done

- **Minimizers:** `min sd`, `min quickmin`.  `MinKokkos`/`MinLineSearchKokkos`
  supply the device machinery, so both are small; they reproduce the CPU
  trajectory step for step.  Registering `min_fire_kokkos` in `Install.sh`, which
  was missing, came along with them.
- **Regions:** `plane`, `prism`, `cylinder`, `cone`, `ellipsoid`.  The concrete-type
  list that `fix wall/region/kk` needs now lives once in `region_kokkos_styles.h`;
  the six fixes that only need `match_all_kokkos` (`addforce`, `aveforce`, `efield`,
  `electron/stopping`, `oneway`, `setforce`) test `dynamic_cast<KokkosBase *>`
  instead of a hardcoded style list, so they accept every ported region
  automatically.
- **Fixes (13):** `flow/gauss`, `damping/cundall`, `viscous/nonlinear`,
  `store/force`, `wall/harmonic/outside` (plus the `wall/harmonic/returned` alias),
  `brownian`, `addtorque/atom`, `settorque/atom`, `nvk`, `baoab`, `wall/piston`,
  `propel/self`, `nve/asphere/noforce`.  `propel/self` errors out on QUAT mode and
  `nve/asphere/noforce` on `atom->superellipsoid_flag`, both because
  `AtomVecEllipsoidKokkos` has no bonus data for them.  `baoab` and `brownian` carry
  `skip_tests: kokkos_*` in their YAML references for the same reason `langevin` and
  `gjf` do: the device RNG stream cannot reproduce the CPU Marsaglia stream.  At
  T = 0, where the noise amplitude vanishes, `baoab/kk` is bit-identical to the CPU
  style including the thermostat energy tally.
- **Pair styles on the `pair_kokkos.h` template (23):**
  - FEP soft-core (8): `lj/cut/soft`, `coul/cut/soft`, `lj/cut/coul/cut/soft`,
    `lj/class2/soft`, `morse/soft`, `coul/long/soft`, `lj/cut/coul/long/soft`,
    `lj/charmm/coul/long/soft`.  No `ncoultablebits` anywhere in `src/FEP/`, so the
    long variants evaluate erfc directly and need no table specialisation.
    `fix adapt/fep` needs no special handling: it calls `Pair::reinit()`, which runs
    `init_one()` for every type pair, and each `/kk` style marks `k_params` modified
    there.
  - CORESHELL `/cs` (8): `coul/long/cs`, `born/coul/long/cs`, `buck/coul/long/cs`,
    `lj/cut/coul/long/cs`, `lj/class2/coul/long/cs`, `coul/wolf/cs`,
    `born/coul/wolf/cs`, `born/coul/dsf/cs`.  The six-term B0..B5 erfc series goes
    with its **own** EWALD_P (9.95473818e-1), not `EwaldConst::EWALD_P`, which pairs
    with the A1..A5 series -- using the wrong one is an O(1) error, not a rounding
    one.  In the four styles that combine LJ with long-range Coulomb the CPU lets
    the EPS_EWALD_SQR correction reach the LJ term through a shared `r2inv`;
    `pair_kokkos.h` calls `compute_fpair` without `factor_coul`, so the `/kk`
    variants apply it to the Coulomb term only.  The difference is bounded by
    EPS_EWALD_SQR/rsq and only shows up for pairs with 0 < special_coul < 1 and
    special_lj > 0 (measured 7.7e-12 relative); the eight affected YAML references
    have their `epsilon` raised to about four times that.
  - Singletons (7): `lj/smooth/linear` (and its `lj/sf` alias), `nm/cut/split`,
    `coul/slater/cut`, `born/coul/dsf`, `lj/expand/sphere`, `lj/relres`,
    `lj/charmmfsw/coul/charmmfsh`.
- **Many-body derivations (4):** `sw/mod`, `tersoff/mod/c`, `gran/hooke`,
  `gran/hertz/history`.  Following the package convention, each derives from its own
  CPU class with the parent KOKKOS implementation copied and edited, as
  `tersoff/mod` and `lj/charmm/coul/charmm/implicit` already do.  `tersoff/mod/c`
  rejects the `shift` keyword like the other KOKKOS tersoff styles, so its shift
  reference now skips the accelerator variants.
- **Bonded table styles (4):** `bond table`, `angle table`, `dihedral table`,
  `dihedral table/cut`.  The data `compute_table()` builds on the host goes into
  device views indexed by table number, and `uf_lookup()` becomes a
  `KOKKOS_INLINE_FUNCTION`.  The bond and angle arrays carry one padding element,
  because their spline branch reads `itable+1` one past the end at the last bin,
  where its weight is exactly zero; the dihedral tables are cyclic and wrap instead.
  `dihedral table` carries a device `minimum_image()` and the box data it needs,
  following `bond_harmonic_restrain/kk`.
- **Computes (7):** `ke`, `com`, `inertia`, `gyration`, `entropy/atom`,
  `centro/atom`, `hexorder/atom`.  `GroupKokkos` gains `gyration_kk()` and
  `gyration_tensor_kk()` next to the existing five reductions.  `compute inertia`
  skips the extended-particle term entirely when no ellipsoid/line/tri/body atom
  style is defined, which is when it contributes nothing, and only then falls back
  to the host loop.  `centro/atom` supports `axes yes`: that path needs only cross
  and normalize, not the eigensolver an earlier version of this file assumed.

### Not done, and why

- **`bond special` is not portable.**  Its `compute()` calls `Pair::single()` on the
  host for every bond, and there is no device `single()` in the KOKKOS pair styles
  to call instead.  Porting it means giving every pair style a device `single()`
  -- a framework item, not a style port.
- **Atom styles** (`sph`, `peri`, `edpd`/`mdpd`/`tdpd`, `electron`, `smd`, `apip`,
  `oxdna`) were deliberately deferred: every atom style the rest of Tier 1 needs
  (`sphere`, `ellipsoid`, `dipole`, `charge`, `full`, `molecular`, `atomic`) is
  already ported, and none of these packages has a ported pair/bond/fix style, so a
  `/kk` atom vec would lift the `AtomKokkos::new_avec` hard error but leave the
  physics on the host.  `AtomKokkos::new_avec` hard-erroring still makes them the
  highest-leverage *gating* items in the scan; each needs new masks in
  `atom_masks.h` (`RHO/DRHO/ESPH/DESPH/VEST`, `EDPD_*`, `CC`,
  `ERADIUS/ERVEL/ERFORCE`), none of which exist today.  Do **not** port
  `atom_style template`: `neighbor_kokkos.cpp:88` errors on molecule templates
  because the KK neighbor build reads special bonds from the per-atom `special`
  list only, so porting the vec without fixing the builder would make neighbor
  lists silently wrong.

### CPU bugs this pass turned up

Two were fixed in place, because a `/kk` variant cannot reproduce them:

- `pair lj/smooth/linear` (and its `/omp` clone, and the `lj/sf` alias) scaled
  `fpair` by `factor_lj` but tallied the pair energy unscaled, so with any
  `special_bonds` setting other than 1.0 the reported van der Waals energy did not
  match the forces, nor what `single()` returns.  `pair_kokkos.h` always applies
  `factor_lj` to `compute_evdwl()`.
- `~PairLJCharmmfswCoulCharmmfsh()` restored the LAMMPS Coulomb conversion constant
  before checking `copymode`, so every functor copy going out of scope reset
  `force->qqr2e` and the next `Force::init()` rebuilt `qqrd2e` from the wrong value.
  This is rule 1 of the KOKKOS instructions, and its
  `lj/charmmfsw/coul/long` sibling already had the guard in the right place.

Two more were left alone, and the `/kk` styles reproduce them:

- `pair morse/soft`'s `compute()` does not subtract `offset[][]` from `evdwl`, while
  its `single()` does.
- `compute entropy/atom`'s first histogram bin has `rbinsq[0] = 0`, so an atom with
  a neighbor inside `3*sigma` gets a division by zero and a NaN.

## Tier 2 -- moderate effort, clear payoff

- **QEQ family (6 fixes + `qeq/comb`).** All build a CSR sparse `H` over the neighbor
  list then run CG with a forward comm per iteration -- exactly what the existing
  1753-line `fix_qeq_reaxff_kokkos.cpp` does.  The shared `FixQEq` base means one
  port of the base machinery covers all six.  **The largest block of genuine,
  unblocked device work in the scan.**  `qeq/dynamic` and `qeq/fire` skip CG.
- `fix indent` -- four per-atom force loops (sphere/cylinder/plane/cone).  Only the
  `compute_equal` scalars need host evaluation + broadcast.  Note it has grown to
  ~1140 lines with four geometry variants; this is no longer a small port.
- `nvt`/`npt`/`nph/asphere` as `FixNHAsphereKokkos : FixNHKokkos<DeviceType>`
  overriding `nve_v`/`nve_x`/`nh_v_temp` for angmom (the nh-sphere lesson), plus a
  new `compute temp/asphere/kk`.
- `pppm/stagger/kk` -- `PPPMStaggerKokkos : PPPMKokkos<DeviceType>`, run the existing
  kernels twice with a shifted `shiftone`.  The most realistic KSPACE target.
- `nb3b/harmonic` (copy `pair_sw_kokkos`, delete `twobody`), then `nb3b/screened`
  rides along free.
- `pair atm` -- clean triple loop over one neighbor list, no tables/comm/map; the
  cleanest *new* many-body kernel available.  Needs `ev_tally3` on device.
- `pair lj/sf/dipole/sf` -- pure cutoff dipole, no kspace dependency, so it is the
  honest DIPOLE target (unlike `lj/cut/dipole/long`, which needs an unported kspace).
- `pair lubricate` -- no bonus data; blockers are host-scalar (query `fix deform`/
  `fix wall` each step, copy into the functor; template on `flaglog`/`flagfld`).
- `pair kolmogorov/crespi/z` and `lebedeva/z` -- no normals at all (normal = z), so
  no ILP machinery.  Small hand kernels; the force is anisotropic, so they do *not*
  fit `compute_fpair`.
- `fix msst`, `fix rattle` (subclass of the ported `FixShake`), `fix wall/table`,
  `wall/ees`, `wall/reflect/stochastic`, `fix gld`, `fix ffl`.
- `compute rdf` (neighbor-list histogram with `atomic_add` into a device bin array),
  `compute sna/atom` (`sna_kokkos_impl.h` already puts the descriptor on device).
- `pair vashishta/table`, `pair gw` (Tersoff-shaped; copy `pair_tersoff_kokkos`).
- `pair coul/exclude` -- needs a hand kernel, not the template: the CPU loop
  deliberately has no `rsq < cutsq` guard and skips `sbmask(j) == 0`.

## Tier 3 -- hard, but the only things that move flagship workloads

- **`pair granular`** (the modern unified granular style, 52 example decks).  Not a
  port, a redesign: physics lives in a runtime-composed hierarchy of virtual
  `GranSubMod` objects invoked per contact through `model->calculate_forces()`.
  Device execution needs the whole sub-model set converted to compile-time tags or
  an integer switch.
- **`pppm/disp`** (8860 lines).  Second grid family plus up to seven dispersion
  bricks for arithmetic mixing plus a ragged `density_brick_none` path.  Mitigations:
  a0..a6 share `gc6`/`fft*_6`, and `poisson_ik`/`poisson_peratom` are already
  parameterized on `FFT3d*`.  **Realistic scoping: `function[0]` (coulomb) +
  `function[1]` (geometric mixing) only, hard-error on arithmetic/none.**  Even that
  is ~2x the existing `pppm_kokkos.cpp`.
- `pair lj/long/coul/long`, `buck/long/coul/long` -- runtime `ewald_order` bitmask
  picks four-plus kernel variants, two table sets, and long-range dispersion changes
  `special_lj` semantics.  Pointless before `pppm/disp/kk`.
- `pair gayberne` then `resquared` -- `AtomVecEllipsoidKokkos` and
  `MathExtraKokkos::quat_to_mat` exist, so the infrastructure is there; the cost is
  per-pair 3x3 algebra and register pressure.  Do gayberne first.
- **INTERLAYER ILP family** -- gated on one piece of machinery: `ILP_neigh()` (a
  private 3-nearest-neighbor list built from the pair list) plus `calc_normal()`
  (per-atom normals and full derivative tensors).  Once `ilp/graphene/hbn/kk` exists,
  `saip/metal`, `ilp/tmd`, `kolmogorov/crespi/full`, `aip/water/2dm` and
  `saip/metal/tmd` ride along nearly free.
- `pair thole` + `lj/cut/thole/long` -- `atom->map(drudeid[i])` + `closest_image` per
  atom; precompute partner indices into a device view once per reneighbor.
- `pair hbond/dreiding/lj` (then `/morse` free) -- walks bond topology via
  `atom->map(klist[kk] + tagprev)`, needs molecule-template `tagprev` offsets and a
  3-body angular force.
- `pair eam/cd`, `meam/spline`, `local/density`, `eim`, `dispersion/d3` -- all need
  the same missing piece, **device pair forward+reverse comm** for a per-atom density
  and its derivative.  Build that helper once and all five become moderate.
- `fix ttm` and `ttm/grid`.
- `ewald/kk` -- deferred as O(N^1.5), but technically the *easiest* kspace here (two
  `parallel_reduce` over kmax x nlocal plus a force kernel), and it has `/gpu`+`/omp`.
  Worth keeping on the list as a cheap device-residency win for small/medium systems.

## Do NOT port

**No device work to offer.** `compute pe` and `compute pressure` (both have
`datamask_read = datamask_modify = EMPTY_MASK`: no per-atom loop at all -- the
scalar is the temperature scalar plus the global `virial[6]`, already reduced
on-device by the force styles) - `fix ave/time`, `fix vector`, `fix print`,
`fix set`, `fix store/state`, `compute property/atom` - all `Command` and `Dump`
styles (setup-time or host I/O).

**Meta-styles over other styles' host pointers** -- a `/kk` would be a fake port:
`compute reduce`, `reduce/region`, `slice`, `global/atom`, `fix pair`, `fix ave/atom`.
These fall out of Tier 0 item 4, not out of individual ports.

**Inherently host: atom count or topology changes.** `fix deposit`, `evaporate`,
`pour`, `append/atoms`, and all of MC (`gcmc`, `widom`, `atom_swap`, `bond_break`,
`bond_create`, `charge/regulation`).  `fix tfmc` additionally has an unbounded
per-atom rejection while-loop.

**Inherently serial algorithms.** `fix thermal/conductivity` and `fix viscosity`
(Mueller-Plathe sorted insertion + `MPI_Allreduce(MAXLOC)`) - `compute born/matrix`,
`stress/mop`, `fabric` (all need the host `pair->single()` virtual, which has no
`pair_kokkos.h` counterpart) - `fix numdiff`, `numdiff/virial` (an `MPI_Allreduce`
per atom per dimension) - `fix gle` (dense matrix product over packed velocities).

**`fix move`, `fix deform/pressure`** -- box-deformation logic and `xoriginal` as a
plain `double **`; host-side scalar work.  **`fix restrain`** -- portable now that a
device map exists, but the restraint list is typically under 100 entries, so there is
no GPU payoff.

**Whole `TALLY` package** -- built on `Pair::add_tally_callback`, a host callback
invoked per-pair inside the CPU pair loop.  `grep -rn "tally_callback" src/KOKKOS/`
returns nothing.  Structurally incompatible.

**The SPIN package -- all of it.**  `atom_vec_spin_kokkos` exists, which makes SPIN
look like a ready target.  It is not.  `fix nve/spin`'s symplectic update is an
explicitly *sequential* sweep: `ComputeInteractionsSpin(i)` recomputes atom i's
neighbor sum, `AdvanceSingleSpin(i)` writes `sp[i]`, forward then backward, with
`comm->forward_comm()` **inside** the loop and a red-black sectoring scheme for MPI,
because each atom must see its neighbors' *updated* spins.  A device port needs a
different algorithm (real graph coloring), not a kernel translation.
`fix langevin/spin` has `setmask()` returning 0 and is called from inside that sweep,
so it cannot be ported independently.  And all six `pair_spin_*` styles are unported,
so even a perfect `precession/spin/kk` would still force a full `sp`/`fm` host
round-trip every step.  All-or-nothing, and the "all" is an algorithm redesign.  The
same reasoning rules out `min spin`, `min spin/cg` and `min spin/lbfgs`.

**External libraries.** `pair lepton` and `bond`/`angle`/`dihedral lepton` (host JIT
expression trees) - `hdnnp` (n2p2) - `quip` (QUIP/GAP) - `voronoi/atom` (voro++) -
`ptm/atom` - `colvars`, `plumed`, `mdi`, `scafacos`, `kim`, `mbx`.

**`kspace msm` and `msm/cg`** -- ragged `double ****` multigrid hierarchy with an
`MPI_Comm` per level, and its niche (non-periodic) is not a GPU workload.  This also
makes `pair coul/msm`, `lj/cut/coul/msm` etc. pointless despite their `/gpu` variants.

**ELECTRODE** (`ewald/electrode`, `pppm/electrode`, `fix electrode/conp`) -- blocked
on the interface, not the kspace: `ElectrodeKSpace::compute_matrix` builds a dense
N_elec x N_elec capacitance matrix for host LAPACK inversion.  Even a perfect
`pppm/electrode/kk` grid path leaves the fix host-bound.

**`atom_style body`** (polymorphic `Body *bptr` plus ragged per-particle payloads) and
**`atom_style template`** -- `neighbor_kokkos.cpp:88` errors on molecule templates
because the KK neighbor build reads special bonds from the per-atom `special` list
only.  Porting the vec without fixing the builder would make neighbor lists
**silently wrong**.

**`fix rigid` (big) and `rigid/nve|nvt|npt|nph`** -- bodies replicated on all procs
with 14 `MPI_Allreduce` over `nbody` per step.  `rigid/small` is the recommended
style and is already fully ported.

**`fix npt/kk` and `nph/kk` are COMPLETE full ports** -- every integrator hook
(`nh_v_press`, `nve_v`, `nve_x`, `nh_v_temp`) and `remap()` (via
`DomainKokkos::x2lamda`/`lamda2x`) runs as a device kernel, including the triclinic
path.  Nothing remains to do on the Nose-Hoover barostat itself.

**`pair airebo`/`airebo/morse`/`airebo/bc`, `rebo`, `comb`, `comb3`** -- variable-depth
bond-path search, torsion, bicubic spline patches, self-consistent QEq.  No
accelerator package in LAMMPS has ever ported AIREBO (OPENMP/INTEL only).

**`pair srp`** -- needs `fix SRP` pseudo bond-particles; `onetwoexclude()` mutates
`firstneigh` in place on the host.

**`fix ilves/omp`** -- per-cluster OpenMP threading was benchmarked and rejected (net
slowdown); the same reasoning argues against a KOKKOS thread-parallel variant.

**r-RESPA** -- four independent blockers: no `RespaKokkos`; every KOKKOS pair style
sets `respa_enable = 0`, so there is nothing to accelerate; `FixRespa::f_level` is a
private `double ***` needing DualView plus exchange packing; and there are no `/kk`
respa neighbor lists.  The largest single item in this assessment -- recommend
against without a concrete user driving it.

## Neighbor: what the raw counts actually mean

Do **not** track `NStencilStyle` or `NPairStyle` as a gap.

- `NStencilStyle` is 16 base styles and 0 `/kk`, and always will be: KOKKOS reuses
  the host-computed stencils (see "Framework facts" above).  A pure counting artifact.
- Of the ~150 non-kk `NPairStyle` registrations: roughly half are `/omp` variants
  (KOKKOS threads via Kokkos); many are `atomonly` micro-variants that
  `NPairKokkos<DeviceType,HALF,NEWTON,GHOST,TRI,SIZE>` folds into one templated
  kernel; every `nsq` variant is deliberately rejected (`neighbor_kokkos.cpp:125` --
  N-squared build is never right on a GPU); `multi/old` was removed from LAMMPS
  entirely; and the respa lists cannot be requested by anything (see r-RESPA above).
- `NBinStyle` is three real base styles (`standard`, `multi`, `ssa`) plus `intel`,
  with two ported (`kk`, `ssa/kk`).

The honest gap list is **one item: `multi` binning** (Tier 0 item 1), plus
`skip/size/off2on` and `skip/size/off2on/oneside` as a distant second (used by
granular hybrid setups and by `fix wall/gran/region` / GRANSURF; `NPairSkipKokkos`
has no OFF2ON or ONESIDE variant).
