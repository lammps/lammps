# KOKKOS Porting Guide: Templates and Lessons

How-to reference for porting LAMMPS pair, fix, and compute styles to the KOKKOS package:
template skeletons plus the accumulated lessons from past porting work, by theme. The
always-on rules (copymode guard, virtual allocate, `AT::` vs `DAT::`, datamask/sync,
full-ports-only, doc checklist) live in `.github/instructions/kokkos.instructions.md` and
are not repeated here; remaining work is in `.github/dev-docs/kokkos-porting-backlog.md`.

## Porting a pair style

Most simple pairwise styles (no three-body terms, no many-body embedding) fit the
`pair_kokkos.h` dispatch template. The KOKKOS class implements:

| Function | COUL_FLAG=0 | COUL_FLAG=1 |
|---|---|---|
| `compute_fpair` | required | required (VdW part) |
| `compute_evdwl` | required | required |
| `compute_fcoul` | return 0 | required |
| `compute_ecoul` | return 0 | required |

The template automatically handles half/half-thread/full neighbor-list dispatch, force
reduction, and energy/virial accumulation. Canonical examples: `pair_born_kokkos` and
`pair_lj_cut_kokkos` in `src/KOKKOS/`; for short-range Coulomb (Wolf/DSF/cut) follow
`pair_lj_cut_coul_dsf_kokkos`; for long-range tables `pair_lj_cut_coul_long_kokkos`
(uses `init_tables()`).

### Header skeleton (COUL_FLAG=0)

```cpp
#ifdef PAIR_CLASS
PairStyle(<name>/kk,     Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/device,Pair<Name>Kokkos<LMPDeviceType>);
PairStyle(<name>/kk/host,  Pair<Name>Kokkos<LMPHostType>);
#else
#ifndef LMP_PAIR_<NAME>_KOKKOS_H
#define LMP_PAIR_<NAME>_KOKKOS_H
#include "pair_kokkos.h"
#include "pair_<name>.h"
#include "neigh_list_kokkos.h"
namespace LAMMPS_NS {
template<class DeviceType>
class Pair<Name>Kokkos : public Pair<Name> {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  Pair<Name>Kokkos(class LAMMPS *);
  ~Pair<Name>Kokkos() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_<name> {
    KOKKOS_INLINE_FUNCTION params_<name>() { /* zero all */ }
    KOKKOS_INLINE_FUNCTION params_<name>(int) { /* zero all */ }
    KK_FLOAT cutsq, /* ... pair-specific params ... */;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(...) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(...) const { return 0; }

  Kokkos::DualView<params_<name>**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_<name>**,...>::t_dev_const_um params;
  params_<name> m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  Kokkos::DualView<KK_FLOAT**,DeviceType> k_cutsq;
  typename Kokkos::DualView<KK_FLOAT**,...>::t_dev d_cutsq;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;
  Kokkos::DualView<KK_FLOAT*,DeviceType> k_eatom, k_vatom;
  int neighflag, newton_pair;
  int nlocal, nall, eflag, vflag;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,FULL,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALF,false>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,true>;
  friend struct PairComputeFunctor<Pair<Name>Kokkos<DeviceType>,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    FULL,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALF,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<Pair<Name>Kokkos<DeviceType>,
    HALFTHREAD,CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<0>>(Pair<Name>Kokkos<DeviceType>*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<Pair<Name>Kokkos<DeviceType>>(
    Pair<Name>Kokkos<DeviceType>*);
};
}
#endif
#endif
```

For `COUL_FLAG=1`, add `compute_fcoul` and `compute_ecoul`, add `special_coul`, `qqrd2e`,
and the `q` array, and change `CoulLongTable<0>` to `CoulLongTable<1>` in the friend
declarations.

### .cpp skeleton

```cpp
template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::Pair<Name>Kokkos(LAMMPS *lmp) : Pair<Name>(lmp)
{
  respa_enable = 0;
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

template<class DeviceType>
Pair<Name>Kokkos<DeviceType>::~Pair<Name>Kokkos()
{
  if (copymode) return;
  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in; vflag = vflag_in;
  if (neighflag == FULL) no_virial_fdotr_compute = 1;
  ev_init(eflag,vflag,0);
  if (eflag_atom) { /* reallocate k_eatom */ }
  if (vflag_atom) { /* reallocate k_vatom */ }
  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);
  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal; nall = atom->nlocal + atom->nghost;
  newton_pair = force->newton_pair;
  special_lj[0..3] = ...;
  NeighListKokkos<DeviceType>* k_list = ...;
  ...
  EV_FLOAT ev = pair_compute<Pair<Name>Kokkos<DeviceType>,
    CoulLongTable<COUL_FLAG>>(this, k_list);
  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) { virial[0..5] += ev.v[0..5]; }
  if (eflag_atom) { k_eatom.template modify<DeviceType>(); k_eatom.sync_host(); }
  if (vflag_atom) { k_vatom.template modify<DeviceType>(); k_vatom.sync_host(); }
  if (vflag_fdotr) pair_virial_fdotr_compute(this);
  copymode = 0;
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::allocate()
{
  Pair<Name>::allocate();
  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  k_params = Kokkos::DualView<params_<name>**,...>("Pair<Name>::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

template<class DeviceType>
void Pair<Name>Kokkos<DeviceType>::init_style()
{
  Pair<Name>::init_style();
  // Request full neighbor list
  neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_GHOST);
}

template<class DeviceType>
double Pair<Name>Kokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = Pair<Name>::init_one(i,j);
  k_params.view_host()(i,j).<field> = static_cast<KK_FLOAT>(<field>[i][j]);
  /* ... fill all params fields ... */
  k_params.view_host()(j,i) = k_params.view_host()(i,j);
  if (i < MAX_TYPES_STACKPARAMS+1 && j < MAX_TYPES_STACKPARAMS+1)
    m_params[i][j] = m_params[j][i] = k_params.view_host()(i,j);
  m_cutsq[j][i] = m_cutsq[i][j] = static_cast<KK_FLOAT>(cutone*cutone);
  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone*cutone;
  k_cutsq.modify_host();
  k_params.modify_host();
  return cutone;
}

// Explicit instantiation at the bottom of the .cpp file
namespace LAMMPS_NS {
template class Pair<Name>Kokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class Pair<Name>Kokkos<LMPHostType>;
#endif
}
```

## Porting a fix style

Key structural differences from pair-style porting:

1. **No dispatch template.** There is no fix equivalent of `pair_kokkos.h`; each fix
   implements its lifecycle hooks directly with `Kokkos::parallel_for` /
   `Kokkos::parallel_reduce`. Primary template for a simple per-atom force loop:
   `fix_viscous_kokkos`.
2. **Port only the hooks that matter.** Identify the hooks the base fix implements
   (`grep "void Fix[A-Z]" src/<PKG>/fix_<name>.cpp`); common targets are `post_force`,
   `initial_integrate`, `final_integrate`, `end_of_step`. Hooks that run rarely
   (e.g. `write_restart`) stay on the CPU via the base-class implementation.

### Fix header skeleton

```cpp
#ifdef FIX_CLASS
// clang-format off
FixStyle(<name>/kk,Fix<Name>Kokkos<LMPDeviceType>);
FixStyle(<name>/kk/device,Fix<Name>Kokkos<LMPDeviceType>);
FixStyle(<name>/kk/host,Fix<Name>Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_<NAME>_KOKKOS_H
#define LMP_FIX_<NAME>_KOKKOS_H

#include "fix_<name>.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFix<Name>{};

template<class DeviceType>
class Fix<Name>Kokkos : public Fix<Name> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  Fix<Name>Kokkos(class LAMMPS *, int, char **);
  ~Fix<Name>Kokkos() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFix<Name>, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3_lr x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
};

}

#endif
#endif
```

### Hook pattern

The constructor sets `kokkosable`, `atomKK`, `execution_space`, and the datamasks as in
the pair-style skeleton above; each ported hook follows the sync/views/kernel/modified sequence:

```cpp
void Fix<Name>Kokkos<DeviceType>::post_force(int /*vflag*/)
{
  atomKK->sync<DeviceType>(datamask_read);
  x    = atomKK->k_x.view<DeviceType>();
  f    = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();
  int nlocal = atom->nlocal;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFix<Name>>(0, nlocal), *this);
  atomKK->modified<DeviceType>(datamask_modify);
}
```

## Lessons: global scalars, reductions, and the host/device split

- **Global scalar parameters** set by `settings()` (not per-pair coefficients) are plain
  `double` base members. Copy them to `KK_FLOAT` members in `compute()` before the
  kernel launch; kernels read them through `this` (the functor copy).
- **Plain int/flag base members are read on device for free** -- captured by value with
  the functor object (like `groupbit`). Only **pointer** members (`boxlo`, `prd`, ...)
  must be copied into device-friendly scalar arrays (`m_boxlo[3]`, ...); host pointers
  cannot be dereferenced on device.
- **Thermostats/barostats: global scalar on the host, per-atom work on the device.**
  Global reductions (`group->ke()`/`vcm()`/`mass()`, a Compute's `compute_scalar()`, a
  `RanMars` draw) stay on the host -- but source the kinetic energy from a *kokkosable*
  temperature compute (`temp/kk`) so the reduction runs on the GPU and only one scalar
  comes back. The per-atom rescale MUST be a device kernel; delegating it to the CPU
  base class is the fake-port trap.
- **Group-average force fixes** (the `aveforce` pattern) need a device
  `parallel_reduce` plus an MPI allreduce before the per-atom apply kernel.

## Lessons: tapering and switching functions

- MDF-style tapers (`buck/mdf` and relatives): `compute_fpair` receives `rsq` and takes
  `r = sqrt(rsq)` itself; inline the taper behind `if (rsq > cut_inner_sq)` and store
  `cut_inner`, `cut_inner_sq`, and the outer `cutsq` in the per-pair `params_*` struct.
- Polynomial tapers (`coul/shield`): inline the taper and its derivative via Horner's
  method with the same hardcoded coefficients as the CPU base, directly inside
  `compute_fcoul`/`compute_ecoul`. No LAMMPS utility function is needed on device.

## Lessons: custom kernels for orientation-dependent styles

Styles that read per-atom quaternions and accumulate torques (`ylz`, `gayberne`, ...) do
not fit the `pair_kokkos.h` template. Write a custom kernel following
`pair_lj_cut_dipole_cut_kokkos`:

- Use `TagPairXXXKernel<NEIGHFLAG,NEWTON_PAIR,EVFLAG,STACKPARAMS>` tags; override
  `compute()` to launch the right instantiation with `Kokkos::parallel_for/reduce`.
- Accumulate torques atomically via a view with `Kokkos::MemoryTraits<AtomicF<NEIGHFLAG>::value>`.
- Access ellipsoid bonus data (quaternions) via `avecKK->k_bonus.view<DeviceType>()`,
  where `avecKK` is an `AtomVecEllipsoidKokkos*` cast from
  `atom->style_match("ellipsoid")`; convert with `MathExtraKokkos::quat_to_mat()`.
- Include `TORQUE_MASK | ELLIPSOID_MASK` in `datamask_read`, `TORQUE_MASK` in `datamask_modify`.
- Store the neighbor-list arrays explicitly (`d_neighbors`, `d_numneigh`, `d_ilist`)
  instead of relying on the `pair_compute` helper.
- **Force sign convention:** the base CPU code accumulates `f[i][0] += delx * fpair`
  with `fpair = F/r` and `delx = xj - xi`; the KOKKOS kernel must match exactly. For
  `ylz` the force on `i` due to `j` is `-dU/dr * rhat` with `rhat` pointing from `i` to
  `j`; the force on `j` is the Newton-3 partner `-force_on_i`.
- **Device math subtlety:** for exponents that are themselves functions of distance
  (`lj/pirani`: `pow(rx, n(rx))`), `Kokkos::pow(base, exp)` is the device-portable form.

## Lessons: per-atom data beyond x/v/f

- **Molecule-id filtering** (`coul/shield` skips same-layer pairs): declare a per-atom
  `molecule` view (`typename AT::t_tagint_1d_randomread`), assign it from
  `atomKK->k_molecule.view<DeviceType>()` in `compute()`, and test
  `if (molecule(i) == molecule(j)) return 0;` in the kernel. (This view was the case
  that motivated the AT-not-DAT rule.)
- **Per-atom radius with a global cutoff ceiling** (`lj/cut/sphere`): `cutsq` stores the
  *maximum possible* squared cutoff per type pair (set in `init_one()`) so the neighbor
  list is conservative; the kernel computes the actual cutoff `(sigma_ij + radius[i] +
  radius[j])^2` on the fly and returns zero force beyond it. Same `typename AT::` rule.

## Lessons: thin-wrapper styles

A global-coefficient variant (`coul/cut/global`) needs no kernels of its own: subclass
the ported parent (`PairCoulCutGlobalKokkos<DeviceType> : PairCoulCutKokkos<DeviceType>`),
override only `coeff()` (argument-count enforcement) and `extract()`, and point the
`PairStyle` macro entries at the new class; the parent handles all kernel logic.

## Lessons: wall styles (FixWall subclasses)

All `FixWall` subclasses share the generic `FixWall::post_force` loop. The KOKKOS port
(`fix_wall_lj93_kokkos` is the template) replaces that loop with a Kokkos functor; each
subclass then needs only `precompute(int)` (per-wall coefficients) and `wall_particle(int
m, int which, double coord)` (single-atom force/energy at a given wall distance).

- **`private:` must become `protected:`** for the base coefficient arrays (`coeff1[]`,
  ..., `offset[]`) that the KK `precompute()` reads.
- **No coefficients, no DualViews:** `wall/harmonic` computes `dr = cutoff - delta`
  inline, so its `precompute()` is empty; the small fixed base arrays (`epsilon[]`,
  `cutoff[]`, 6 elements) are read as scalar casts in the functor. Only `k_vatom` keeps
  the DualView pattern.
- **`powint()` replacement** (`wall/lj1043`): store the needed coefficient array in its
  own DualView and expand integer powers by hand in the functor:

  ```cpp
  KK_FLOAT dc4inv = 1.0 / (delta + d_coeff4(m));
  KK_FLOAT dc4inv3 = dc4inv * dc4inv * dc4inv;      // replaces powint(..., -3)
  KK_FLOAT dc4inv4 = dc4inv3 * dc4inv;              // replaces powint(..., -4)
  ```

- `exp()` becomes `Kokkos::exp()` (`wall/morse`); copy the protected base parameters
  (`alpha[]`, `sigma[]`, `epsilon[]`) into DualViews in `precompute`.
- The KK `post_force` override only (1) reallocates `k_vatom` if `vflag_atom`, (2) calls
  the base `post_force(vflag)` -- which virtual-dispatches into the KK `wall_particle`
  -- and (3) syncs `k_vatom` back to host. Virial accumulation in `wall_particle` uses
  the `result[7..12]` reduction slots, flushed to `virial[0..5]` after `parallel_reduce`.

## Lessons: integrator hierarchies (the NH-sphere pattern)

- **Inherit from the KOKKOS mid-class, not the CPU mid-class.** The instinct is
  `FixNHSphereKokkos : FixNHSphere` (as the wall styles do), but `FixNHSphere::nve_v()`
  calls the CPU `FixNH::nve_v()` internally, losing all Kokkos parallelism. Correct:
  `FixNHSphereKokkos<DeviceType> : FixNHKokkos<DeviceType>` (brings in the device
  `nve_v`, `nve_x`, `nh_v_temp`, `nh_v_press`); override just the sphere-specific
  methods to call the base version then launch the omega-update kernel, and duplicate
  the `FixNHSphere` constructor logic (omega/radius checks, `inertia`, `disc`).
- **Name collisions:** `FixNH` has `protected: double omega[6]` (barostat DOFs); name
  the per-atom angular-velocity view `omega_kk`.
- **Base views stay valid across the virtual-call chain:** after `FixNHKokkos::nve_v()`
  returns, `this->mask` / `this->rmass` still hold the current device views (reference-
  counted handles); the omega kernel reuses them, syncing only OMEGA/TORQUE/RADIUS masks.
- **Thin wrappers must create the kk helper:** `nvt/sphere/kk` and siblings create
  `temp/sphere/kk` -- not `temp/sphere` (correct but forces a host sync every thermostat
  step) and not `temp/kk` (wrong temperature: ignores angular velocity).
- **DLM dipole integrator unsupported on device** (matrix-product chain, rarely used):
  support only `dlm_flag == 0` and raise an error for `update dipole dlm`.
- `compute temp/sphere`: sphere particles always provide `rmass`, so the `RMASS`
  template parameter of `compute temp` is unnecessary -- use a `MODE` parameter
  (ROTATE-only vs ALL) instead. `dof_compute()` and `init()` stay CPU-side (infrequent).

## Lessons: temperature-compute port pattern

Template for porting any `compute temp/*` (mirrors `compute_temp_kokkos`): a
`parallel_reduce` into the 6-component `s_CTEMP` struct, templated on `RMASS`, the
style's flags folded into the per-atom terms.
- **A biased temperature (`tempbias = 1`) MUST implement the device bias methods.**
  `Compute::remove_bias_all_kk()` / `restore_bias_all_kk()` are empty stubs in the base;
  with `kokkosable = 1` and no override, a KK thermostat silently skips bias removal --
  wrong results, no error. Implement both as device kernels, plus host-facing
  `remove_bias_all` / `restore_bias_all` that call the `_kk` version then
  `atomKK->sync(Host, V_MASK)`, and override `reapply_bias_all` (used by `velocity ...
  bias yes`) when the CPU base overrides it. Pattern: `compute_temp_deform_kokkos`.
- Store the bias in a device view (`typename AT::t_kkfloat_1d_3 vbiasall`, shadowing the
  base `double **` member, which stays null); `maxbias = 0` in the ctor, regrow when
  `atom->nmax > maxbias`.
- The base destructor needs `if (copymode) return;` -- the compute is copied as a Kokkos
  functor and the copy's destruction chains into the base destructor.
- `dof_compute()` stays on the host (cheap, infrequent).

Binned/profile temperatures (per-bin streaming velocity subtraction) additionally:

- Flip the base `private:` members the KK class needs to `protected:`.
- Keep the per-bin reduction's MPI round-trip ON THE HOST, by design: a scatter kernel
  does `atomic_add` into a device accumulator, copy out to the base array, run the
  *identical* host `MPI_Allreduce` + divide, and copy the averages back to a device
  view. The host reduce is what keeps CPU/KK bit-identical; the bin count is small, so
  the round-trip is cheap, once per `compute_scalar`.
- Assign per-atom bins on device for orthogonal boxes; fall back to the host
  `bin_assign()` (needs lambda-space coordinates) for triclinic. Diagnostic outputs
  (`compute_array`) may sync to host and delegate to the base.

Validation: run a thermostatted MD driven by the compute and require bit-identical CPU
vs `-sf kk` trajectories (exercises the scalar, bias remove/restore, and vector paths).

## Lessons: per-atom-output computes

Template for a per-atom vector-output compute whose result comes purely from per-atom
data (`ke/atom`; pattern: `compute_coord_atom_kokkos`):

- Keep the base `double *` output pointer (make the base member `protected`) and add
  `DAT::ttransform_kkfloat_1d k_<out>` plus `typename AT::t_kkfloat_1d d_<out>`; the
  `ttransform_` view auto-converts host double <-> device `KK_FLOAT`, so it also works
  in single/mixed-precision builds. Grow with `memoryKK->destroy_kokkos/create_kokkos
  (k_<out>, <out>, nmax, ...)`; set `vector_atom` (or `array_atom`).
- After the kernel: `k_<out>.modify<DeviceType>(); k_<out>.sync_host();` so host
  consumers (dumps, `compute reduce`, ...) see current data. `copymode` guards in both
  destructors; `destroy_kokkos` nulls the pointer, making the base `memory->destroy` a
  safe no-op. Scalars the kernel needs (e.g. `mvv2e`) are plain members set in
  `compute_peratom()`.

## Lessons: barostat-fix port pattern

`press/berendsen/kk` is the cleanest template for a "global scalar on host, per-atom
work on device" barostat/thermostat fix:

- **The internal `temp` compute is auto-promoted to `temp/kk` by the suffix machinery**
  (`Modify::add_compute` defaults to `trysuffix = 1`), so under `-sf kk` the KK fix does
  not recreate it. The fake port's one real flaw was NOT ensuring a kk temperature (a
  non-kk `temp` reads `atom->v` on the host every `end_of_step`); the explicit-suffix
  hole without `-sf kk` is closed by the A+C mitigation in the instructions file.
- **`compute pressure` stays non-kk and that is correct** (see the do-NOT-port list in
  the backlog: `EMPTY_MASK`, no per-atom data).
- **`remap()` becomes device-resident for free via virtual dispatch:** `domain` is a
  `DomainKokkos` instance under the KOKKOS package, so the base `domain->x2lamda(nlocal)`
  / `lamda2x(nlocal)` already run on device (each syncs `X` itself). Override `remap()`
  only to swap the `dilate partial` host loop for the device group variant
  `x2lamda(nlocal, groupbit)`; make the base `remap()` virtual (a one-line base change).
- **`end_of_step()` mirrors the base but guards the temperature call:**
  `if (temperature->kokkosable) temperature->compute_scalar();` else the explicit
  sync/modified dance (handles `fix_modify temp` with a non-kk compute). When the fix
  owns no per-atom kernel, `datamask_read = datamask_modify = EMPTY_MASK`, no destructor.

`press/langevin/kk` adds the stochastic-piston and triclinic-flip variations:

- **No temperature compute at all:** the Gronbech-Jensen/Farago barostat uses `NkT/V`
  for the kinetic pressure term and creates its pressure compute as `pressure NULL
  virial`. Do NOT add `F_MASK` to `datamask_read`; the fix reads no per-atom forces.
- **The 6 piston DOFs and random forces stay on the host** (at most 6 Gaussians on rank
  0 + `MPI_Bcast`); `initial_integrate`, `post_force`, and `end_of_step` touch only
  scalars plus the global pressure tensor and are NOT overridden -- only `remap()` and
  `pre_exchange()` are.
- **Copy the triclinic tilt logic verbatim from the base** (including its
  `p_flag[3]->xy` / `[5]->yz` index mapping) so CPU and KK stay identical -- do not
  "fix" base-class quirks inside a port.
- **`pre_exchange()` (box flips) mirrors `fix_nh_kokkos`:** `domainKK->image_flip()` +
  `remap_all()` + `x2lamda()` on device, `atomKK->sync(Host, ALL_MASK)` only around the
  intrinsically host-side `irregular->migrate_atoms()`, then `lamda2x()`. Register the
  hook in `setmask()` only when tilt flips can occur; it fires rarely. Cache `domainKK`
  in the constructor.

## Lessons: random numbers on device

Device RNG uses `Kokkos::Random_XorShift64_Pool<DeviceType>`, with `RandPoolWrap` /
`RandWrap` (`src/KOKKOS/rand_pool_wrap_kokkos.h`, wraps `RanMars` in the pool API) as a
deterministic debug substitute under `#ifndef LMP_KOKKOS_DEBUG_RNG` / `#else` guards
(mirror `gjf`'s constructor/destructor guards):

```cpp
rand_type rand_gen = rand_pool.get_state();
KK_FLOAT r = rand_gen.normal();   // Gaussian deviate
rand_pool.free_state(rand_gen);
```

Choose the strategy by how many draws are needed:

- **A single global factor -> host RNG, bit-identical** (`temp/csvr`): draw the rescale
  factor on rank 0 with the base-class code (`resamplekin`, a few `RanMars` calls),
  `MPI_Bcast`, then apply `v *= lamda` in a device kernel. CPU vs `-sf kk` stays
  bit-for-bit; BIAS works via the device bias methods. No device pool needed.
- **Per-atom draws -> device pool, statistical validation only** (`temp/csld`, `gjf`):
  each atom draws Gaussians on device; the trajectory cannot match the CPU RNG stream,
  so validate that the temperature converges to the target instead. Gotchas: if the base
  constructor consumed the seed into a local variable, re-parse it in the KK
  member-initializer list (`utils::inumeric(FLERR, arg[N], false, lmp) + comm->me`). A
  per-step scratch view (csld's `vhold`) needs no `KokkosBase`/`pack_exchange`; per-atom
  state that migrates with atoms (gjf's half-step velocity `lv`) does.
- **Bias paths that mix saved and fresh velocities** may not map onto
  `remove_bias_all_kk` over `atom->v`; if so, `error->all` on the BIAS option rather
  than silently computing the wrong thing (csld errors out; csvr supports BIAS).

## Lessons: per-atom virial/energy aggregation needs framework work first

`compute stress/atom`, `pe/atom`, and `centroid/stress/atom` are NOT clean full ports and
were deliberately left unported:

- They sum per-atom `vatom`/`eatom` from all seven force-style hierarchies, whose base
  classes expose only host `double **`; the device dual views (`k_vatom`/`d_vatom`) live
  inside each KOKKOS subclass with no uniform accessor.
- They reverse-comm the per-atom result (ghost to owner); no KOKKOS compute does device
  reverse-comm today (needs `KokkosBase` + `pack/unpack_reverse_comm_kokkos` +
  `CommKokkos` support).
- Payoff is marginal anyway: KOKKOS force styles call `k_vatom.sync_host()`
  unconditionally when per-atom virial is requested, so the data reaches the host every
  step regardless.

If revisited, do the framework work first (uniform device per-atom-virial accessor,
device compute reverse-comm, lazy `vatom` host-sync); do NOT ship a host-aggregation
wrapper (a fake port).

## Known limitation: 48 KB team-scratch cap (pair pace/kk)

`pace/kk` and `pace/extrapolation/kk` abort on GPU at high coordination numbers. The
ComputeNeigh kernel requests team scratch of
`team_size (32 on GPU) * maxneigh * sizeof(int)`; Kokkos caps level-0 team scratch at
`sharedMemPerBlock` -- the 48 KB default on every NVIDIA architecture (Kokkos never opts
into the larger per-block dynamic shared memory), so the ceiling (roughly maxneigh 381)
is architecture-independent: a newer NVIDIA GPU does not help. HIP checks the 64 KB LDS
on RDNA (roughly maxneigh 509), so AMD tolerates somewhat more. Symptoms: CUDA
"Requested too much scratch memory on level 0"; HIP "could not find a valid team size".
Unsolved; tracked as issue #5063. Gotchas that defeated two fix attempts:

- `TeamPolicy::team_size_max()` THROWS "could not find a valid team size" instead of
  returning 0 when nothing fits -- you cannot defensively probe an over-budget PerTeam
  policy; the first probe aborts.
- `TeamPolicy::scratch_size_max(0)` does not subtract the kernel's own static shared
  memory, so shrinking team_size to that byte budget still overshoots.
- Candidate designs: request scratch as `PerThread(maxneigh*sizeof(int))` (one thread
  always fits, so team_size_max cannot throw), or drop the shared `inside[]` cache for
  a global `natoms x maxneigh` View, removing the team-size coupling entirely.

Separate GPU restriction (do not conflate with the scratch cap): the `recursive`
evaluator keyword has no device port -- `pair pace/kk` on the GPU requires the `product`
keyword (a controlled `error->all`, independent of maxneigh).

## Debugging GPU-only test failures

Unit-test failures that appear only on a HIP/CUDA build (host backend != device backend)
but pass on Serial/OpenMP-only builds are almost always device/host data-movement bugs,
not physics; the kk math usually matches the CPU base byte-for-byte. Method: reproduce
standalone with CPU `lmp` (no suffix) vs `lmp -k on g 1 -sf kk`, dump `id x y z fx fy fz`
with `dump_modify ... sort id format float %20.12g`, and trace per step (`dump every 1`
over a short run). A step-0 match with growth over steps localizes the bug to the run
loop / integrator data movement. Mirror the force-style test driver exactly: it does
`init_lammps` (`run 0`) THEN `run_lammps` (`fix nve; run 4`), and uses `pair_modify
compute no` + `boundary p p f` (slab) to isolate kspace -- a single run or one-shot will
not reproduce a two-run state bug.

Three GPU-only data-movement bug classes seen (and fixed) in this codebase:

- A `/kk/host` variant calling `commKK->forward_comm_device<DeviceType>(...)` directly:
  the `LMPHostType` instantiation becomes an undefined symbol in `liblammps.so` (breaks
  dlopen for Python/plugins). Fix: explicit host instantiation in `comm_kokkos.cpp`
  under `#ifdef LMP_KOKKOS_GPU`.
- A device-execution-space fix whose host-delegating callbacks do
  `atomKK->modified(Host, mask)` collides with ModifyKokkos's outer
  `modified(Device, mask)`: DualView "concurrent modification" abort. Fix: reconcile
  with `atomKK->sync(execution_space, mask)` after each host modify.
- VerletKokkos's GPU-aware overlap force merge zeroed the host force buffer only when
  `pair_compute_flag` was set; with `pair_modify compute no` + host bonded styles +
  device kspace, stale host force was re-added every step (linear force blow-up). Fix:
  zero it unless the pair style itself runs on host.

Harness tips: `TEST_ARGS=-v` is the first move on any silent-abort or "Only one stdout
capturer can exist at a time" failure -- the drivers echo each captured section as it
completes, so the last echoed line pinpoints where the throw escaped. A style with no
`/kk` variant that cannot run under VerletKokkos needs `skip_tests: kokkos_gpu` in YAML.

Mixed precision is a separate validation target per backend: a style that passes the
mixed-precision fixtures in a HOST-only KOKKOS build can still produce NaNs on a GPU
build when intermediate float-range limits are hit (class example: `exp()` of a large
positive argument overflows in `float` -- reaxff/kk).  Run the `kokkos_gpu_mixed`
fixture on a real device before declaring mixed-precision support clean; if a style is
numerically unfit for single/mixed, add the corresponding `skip_tests` entries instead.

## KOKKOS test build

```bash
cmake -S cmake -B build -C cmake/presets/gcc.cmake -C cmake/presets/most.cmake \
  -D PKG_KOKKOS=on -D Kokkos_ENABLE_OPENMP=on \
  -D DOWNLOAD_POTENTIALS=off -D ENABLE_TESTING=on -G Ninja
cmake --build build -j 4
cd build && ctest -V -R pair
```
