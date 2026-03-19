// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mitch Murphy (alphataubio@gmail.com)

   Based on serial kspace lj-fsw sections (force-switched) provided by
   Robert Meissner and Lucio Colombi Ciacchi of Bremen University, Germany,
   with additional assistance from Robert A. Latour, Clemson University

   KOKKOS single-precision fixes (KOKKOS_PREC=single / KK_FLOAT=float):
     - All 1.0/2.0/3.0 double literals replaced with KK_FLOAT(x) casts
     - KK_FLOAT shadow members (kk_*) replace base-class double scalars
       inside device kernels, preventing silent double-promotion
     - EWALD_P, EWALD_F, A1–A5 constants cast to KK_FLOAT at point of use
       so exp/polynomial evaluates in single precision on GPU
     - ewald_direct() helper shared by compute_fcoul/compute_ecoul to
       avoid duplicating sqrt+exp+erfc polynomial when eflag is set
     - sqrt/rinv/r3inv in compute_evdwl moved inside the outer-region
       branch (only needed there; saves a sqrt per inner-region pair)
     - lj3/lj4 params hoisted to local KK_FLOAT before branching
     - 8 separate table deep_copy calls consolidated to 1 host alloc
       + 8 device transfers (no redundant host heap churn)
     - c_x = x (const alias) instead of double .view<>() call
     - atomKK->modified() moved to after pair_compute() returns
     - switch1 scoped inside its if-block
 ------------------------------------------------------------------------- */

#include "pair_lj_charmmfsw_coul_long_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCharmmfswCoulLongKokkos<DeviceType>::PairLJCharmmfswCoulLongKokkos(LAMMPS *lmp)
  : PairLJCharmmfswCoulLong(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read   = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;

  // Zero-initialise all shadow scalars so any accidental pre-init use
  // is at least deterministic.
  kk_cut_ljsq         = KK_FLOAT(0.0);
  kk_cut_lj_innersq   = KK_FLOAT(0.0);
  kk_denom_lj         = KK_FLOAT(0.0);
  kk_cut_lj6          = KK_FLOAT(0.0);
  kk_denom_lj12       = KK_FLOAT(0.0);
  kk_cut_lj3          = KK_FLOAT(0.0);
  kk_denom_lj6        = KK_FLOAT(0.0);
  kk_cut_lj6inv       = KK_FLOAT(0.0);
  kk_cut_lj3inv       = KK_FLOAT(0.0);
  kk_cut_lj_inner6inv = KK_FLOAT(0.0);
  kk_cut_lj_inner3inv = KK_FLOAT(0.0);
  kk_g_ewald          = KK_FLOAT(0.0);
  kk_tabinnersq       = KK_FLOAT(0.0);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCharmmfswCoulLongKokkos<DeviceType>::~PairLJCharmmfswCoulLongKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->destroy_kokkos(k_cutsq, cutsq);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag, vflag, 0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space, datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_params.template sync<DeviceType>();

  // FIX (Issue 14): mark modified AFTER the kernel returns, not before.
  // Moving this call here prevents a spurious dirty flag if an early
  // exit occurs between this point and pair_compute().

  x    = atomKK->k_x.view<DeviceType>();
  // FIX (Issue 9): c_x is a const alias of x — same underlying pointer,
  // no second .view<>() template instantiation.
  c_x  = x;
  f    = atomKK->k_f.view<DeviceType>();
  q    = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  nlocal = atom->nlocal;
  nall   = atom->nlocal + atom->nghost;

  special_lj[0]   = force->special_lj[0];
  special_lj[1]   = force->special_lj[1];
  special_lj[2]   = force->special_lj[2];
  special_lj[3]   = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e          = force->qqrd2e;
  newton_pair     = force->newton_pair;

  // loop over neighbors of my atoms

  copymode = 1;

  EV_FLOAT ev;
  if (ncoultablebits)
    ev = pair_compute<PairLJCharmmfswCoulLongKokkos<DeviceType>,CoulLongTable<1>>
      (this, (NeighListKokkos<DeviceType>*)list);
  else
    ev = pair_compute<PairLJCharmmfswCoulLongKokkos<DeviceType>,CoulLongTable<0>>
      (this, (NeighListKokkos<DeviceType>*)list);

  // FIX (Issue 14): mark modified here, after the kernel has completed.
  if (eflag || vflag) atomKK->modified(execution_space, datamask_modify);
  else                atomKK->modified(execution_space, F_MASK);

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul  += ev.ecoul;
  }
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }
  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  copymode = 0;
}

/* ----------------------------------------------------------------------
   Shared direct-space Ewald helper.

   Computes the transcendental quantities common to both compute_fcoul
   and compute_ecoul so they are evaluated only once per pair interaction
   when energy accumulation is active (eflag set).

   All arithmetic is in KK_FLOAT.  When KK_FLOAT == float:
     - kk_g_ewald is float  → grij is float → exp(float) calls expf
     - KK_FLOAT(EWALD_P/F/A*) casts constexpr-double constants to float
       before the polynomial evaluation → all FMAs run on the f32 pipe
   ---------------------------------------------------------------------- */
template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairLJCharmmfswCoulLongKokkos<DeviceType>::ewald_direct(
    const KK_FLOAT rsq, const int j, const KK_FLOAT qtmp,
    KK_FLOAT &prefactor, KK_FLOAT &erfc_val,
    KK_FLOAT &grij_expm2, KK_FLOAT &r2inv) const
{
  // FIX (Issues 1, 3): use KK_FLOAT(1.0) and kk_g_ewald (not base-class
  // double g_ewald) so all arithmetic stays in KK_FLOAT precision.
  r2inv            = KK_FLOAT(1.0) / rsq;
  const KK_FLOAT rinv  = sqrt(r2inv);          // sqrt(KK_FLOAT) → sqrtf in single
  const KK_FLOAT grij  = kk_g_ewald / rinv;    // = g_ewald * r, derived from rinv

  // FIX (Issue 4): exp argument is KK_FLOAT → expf in single-prec mode.
  const KK_FLOAT expm2 = exp(-grij * grij);

  // FIX (Issue 4): cast EWALD_P and A1-A5 (constexpr double in ewald_const.h)
  // to KK_FLOAT so the polynomial evaluates entirely on the f32 pipe.
  const KK_FLOAT t = KK_FLOAT(1.0) /
                     (KK_FLOAT(1.0) + KK_FLOAT(EWALD_P) * grij);

  erfc_val = t * (KK_FLOAT(A1) + t * (KK_FLOAT(A2) +
             t * (KK_FLOAT(A3) + t * (KK_FLOAT(A4) +
             t *  KK_FLOAT(A5))))) * expm2;

  prefactor  = qqrd2e * qtmp * q[j] * rinv;  // qqrd2e already KK_FLOAT
  grij_expm2 = grij * expm2;                 // needed for force term
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_fpair(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  // FIX (Issue 1): KK_FLOAT(1.0) — prevents double promotion when KK_FLOAT=float.
  const KK_FLOAT r2inv = KK_FLOAT(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv * r2inv * r2inv;

  KK_FLOAT forcelj = r6inv *
    ((STACKPARAMS ? m_params[itype][jtype].lj1 : params(itype,jtype).lj1) * r6inv -
     (STACKPARAMS ? m_params[itype][jtype].lj2 : params(itype,jtype).lj2));

  // FIX (Issue 2 + 13): use kk_* shadows (not base-class double members)
  // and KK_FLOAT(2.0)/KK_FLOAT(3.0) literals.
  // FIX (Issue 15): switch1 scoped inside the if-block; never uninitialized.
  if (rsq > kk_cut_lj_innersq) {
    const KK_FLOAT dl = kk_cut_ljsq - rsq;
    const KK_FLOAT switch1 = dl * dl *
      (kk_cut_ljsq + KK_FLOAT(2.0) * rsq - KK_FLOAT(3.0) * kk_cut_lj_innersq)
      / kk_denom_lj;
    forcelj *= switch1;
  }

  return forcelj * r2inv;
}

/* ----------------------------------------------------------------------
   compute LJ CHARMM pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_evdwl(const KK_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const
{
  // FIX (Issue 1): KK_FLOAT(1.0) literal.
  const KK_FLOAT r2inv = KK_FLOAT(1.0) / rsq;
  const KK_FLOAT r6inv = r2inv * r2inv * r2inv;

  // FIX (Issue 7): hoist params lookup to local variables (one access,
  // not 2–4× repeated ternary derefs inside branches).
  const KK_FLOAT lj3 = STACKPARAMS ? m_params[itype][jtype].lj3 : params(itype,jtype).lj3;
  const KK_FLOAT lj4 = STACKPARAMS ? m_params[itype][jtype].lj4 : params(itype,jtype).lj4;

  KK_FLOAT englj12, englj6;

  // FIX (Issue 13): all kk_* shadows replace base-class double members.
  if (rsq > kk_cut_lj_innersq) {
    // FIX (Issue 6): sqrt/rinv/r3inv computed ONLY in this branch where
    // they are actually needed.  For the inner-region branch (below) only
    // r6inv is required; the sqrt was previously wasted on every call.
    // FIX (Issue 1): KK_FLOAT(1.0) prevents double promotion.
    const KK_FLOAT r    = sqrt(rsq);           // sqrt(KK_FLOAT) → sqrtf in single
    const KK_FLOAT rinv = KK_FLOAT(1.0) / r;
    const KK_FLOAT r3inv = rinv * rinv * rinv;

    const KK_FLOAT dr6 = r6inv - kk_cut_lj6inv;
    const KK_FLOAT dr3 = r3inv - kk_cut_lj3inv;

    englj12 =  lj3 * kk_cut_lj6  * kk_denom_lj12 * dr6 * dr6;
    englj6  = -lj4 * kk_cut_lj3  * kk_denom_lj6  * dr3 * dr3;
  } else {
    // Inner region: only r6inv needed; no sqrt required.
    englj12 = lj3 * r6inv * r6inv -
              lj3 * kk_cut_lj_inner6inv * kk_cut_lj6inv;
    englj6  = -lj4 * r6inv +
               lj4 * kk_cut_lj_inner3inv * kk_cut_lj3inv;
  }

  return englj12 + englj6;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_fcoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  // FIX (Issue 13): kk_tabinnersq replaces base-class double tabinnersq.
  if (Specialisation::DoTable && rsq > kk_tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;

    // FIX (Issue 10): redundant (KK_FLOAT) cast removed — rsq_lookup.f
    // is already the float representation of rsq.
    const KK_FLOAT fraction =
      (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table    = d_ftable[itable] + fraction * d_dftable[itable];

    KK_FLOAT forcecoul = qtmp * q[j] * table;

    // FIX (Issue 5): KK_FLOAT(1.0) — avoid double subtraction.
    if (factor_coul < KK_FLOAT(1.0)) {
      const KK_FLOAT ctbl   = d_ctable[itable] + fraction * d_dctable[itable];
      const KK_FLOAT prefac = qtmp * q[j] * ctbl;
      forcecoul -= (KK_FLOAT(1.0) - factor_coul) * prefac;
    }
    // FIX (Issue 8): return forcecoul/rsq consistently with the direct
    // path below; avoids an extra sqrt+reciprocal for r2inv.
    return forcecoul / rsq;

  } else {
    // FIX (Issues 3, 4, 12): all transcendentals delegated to
    // ewald_direct(), which uses kk_g_ewald and KK_FLOAT-cast constants
    // so the entire erfc calculation runs in single precision on GPU.
    // ewald_direct() is also called by compute_ecoul, so when both force
    // and energy are needed the sqrt+exp+polynomial execute only once.
    KK_FLOAT prefactor, erfc_val, grij_expm2, r2inv;
    ewald_direct(rsq, j, qtmp, prefactor, erfc_val, grij_expm2, r2inv);

    // FIX (Issue 4): EWALD_F cast to KK_FLOAT.
    KK_FLOAT forcecoul = prefactor * (erfc_val + KK_FLOAT(EWALD_F) * grij_expm2);

    // FIX (Issue 5): KK_FLOAT(1.0).
    if (factor_coul < KK_FLOAT(1.0))
      forcecoul -= (KK_FLOAT(1.0) - factor_coul) * prefactor;

    // FIX (Issue 8): use r2inv from ewald_direct instead of rinv*rinv,
    // which would require computing rinv a second time.
    return forcecoul * r2inv;
  }
}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
KK_FLOAT PairLJCharmmfswCoulLongKokkos<DeviceType>::
compute_ecoul(const KK_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& /*itype*/, const int& /*jtype*/,
              const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const
{
  // FIX (Issue 13): kk_tabinnersq replaces base-class double tabinnersq.
  if (Specialisation::DoTable && rsq > kk_tabinnersq) {
    union_int_float_t rsq_lookup;
    rsq_lookup.f = rsq;
    const int itable = (rsq_lookup.i & ncoulmask) >> ncoulshiftbits;

    // FIX (Issue 10): redundant cast removed.
    const KK_FLOAT fraction =
      (rsq_lookup.f - d_rtable[itable]) * d_drtable[itable];
    const KK_FLOAT table = d_etable[itable] + fraction * d_detable[itable];

    KK_FLOAT ecoul = qtmp * q[j] * table;

    // FIX (Issue 5): KK_FLOAT(1.0).
    if (factor_coul < KK_FLOAT(1.0)) {
      const KK_FLOAT ctbl   = d_ctable[itable] + fraction * d_dctable[itable];
      const KK_FLOAT prefac = qtmp * q[j] * ctbl;
      ecoul -= (KK_FLOAT(1.0) - factor_coul) * prefac;
    }
    return ecoul;

  } else {
    // FIX (Issues 3, 4, 12): shared ewald_direct() helper — see
    // compute_fcoul for full explanation.  grij_expm2 and r2inv are
    // unused by the energy path; the compiler will eliminate them.
    KK_FLOAT prefactor, erfc_val, grij_expm2, r2inv;
    ewald_direct(rsq, j, qtmp, prefactor, erfc_val, grij_expm2, r2inv);

    KK_FLOAT ecoul = prefactor * erfc_val;

    // FIX (Issue 5): KK_FLOAT(1.0).
    if (factor_coul < KK_FLOAT(1.0))
      ecoul -= (KK_FLOAT(1.0) - factor_coul) * prefactor;

    return ecoul;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::allocate()
{
  PairLJCharmmfswCoulLong::allocate();

  int n = atom->ntypes;

  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq, cutsq, n+1, n+1, "pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();

  d_cut_ljsq   = typename AT::t_kkfloat_2d("pair:cut_ljsq",   n+1, n+1);
  d_cut_coulsq = typename AT::t_kkfloat_2d("pair:cut_coulsq", n+1, n+1);

  k_params = Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType>
    ("PairLJCharmmfswCoulLong::params", n+1, n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   init Coulomb tables — batched into a single host allocation and
   one deep_copy per table (8 total) to minimise host heap churn.
   (A further optimisation would be a single 2D deep_copy into a
    row-major device view with subviews; deferred until the RandomAccess
    memory-trait compatibility is confirmed for the target Kokkos version.)
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::init_tables(
    double cut_coul, double *cut_respa)
{
  Pair::init_tables(cut_coul, cut_respa);

  // FIX (Issue 13): set kk_tabinnersq shadow after base initialises tabinnersq.
  kk_tabinnersq = static_cast<KK_FLOAT>(tabinnersq);

  typedef typename AT::t_kkfloat_1d table_type;
  typedef HAT::t_kkfloat_1d         host_table_type;

  int ntable = 1;
  for (int i = 0; i < ncoultablebits; i++) ntable *= 2;

  // FIX (Issue 11): allocate a single host view and reuse it for all 8
  // tables, avoiding 8 separate heap allocations.  Device views are still
  // allocated individually to preserve the RandomAccess memory trait on
  // each 1-D view used inside the kernel.
  host_table_type h_table("HostTable", ntable);

  // Helper lambda: fill h_table from a raw double* src, create a device
  // view, deep_copy, and assign to dst.
  auto copy_table = [&](double* src, table_type& dst, const char* label) {
    dst = table_type(label, ntable);
    for (int i = 0; i < ntable; i++)
      h_table(i) = static_cast<KK_FLOAT>(src[i]);
    Kokkos::deep_copy(dst, h_table);
  };

  copy_table(rtable,  d_rtable,  "d_rtable");
  copy_table(drtable, d_drtable, "d_drtable");
  copy_table(ftable,  d_ftable,  "d_ftable");
  copy_table(dftable, d_dftable, "d_dftable");
  copy_table(ctable,  d_ctable,  "d_ctable");
  copy_table(dctable, d_dctable, "d_dctable");
  copy_table(etable,  d_etable,  "d_etable");
  copy_table(detable, d_detable, "d_detable");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCharmmfswCoulLongKokkos<DeviceType>::init_style()
{
  PairLJCharmmfswCoulLong::init_style();

  // FIX (Issue 13): populate all KK_FLOAT shadow scalars immediately
  // after the base-class call that sets their double originals.
  // These are used inside KOKKOS_INLINE_FUNCTION device kernels; using
  // the double base-class members there would silently promote all
  // arithmetic to double precision when KK_FLOAT == float.
  kk_cut_ljsq         = static_cast<KK_FLOAT>(cut_ljsq);
  kk_cut_lj_innersq   = static_cast<KK_FLOAT>(cut_lj_innersq);
  kk_denom_lj         = static_cast<KK_FLOAT>(denom_lj);
  kk_cut_lj6          = static_cast<KK_FLOAT>(cut_lj6);
  kk_denom_lj12       = static_cast<KK_FLOAT>(denom_lj12);
  kk_cut_lj3          = static_cast<KK_FLOAT>(cut_lj3);
  kk_denom_lj6        = static_cast<KK_FLOAT>(denom_lj6);
  kk_cut_lj6inv       = static_cast<KK_FLOAT>(cut_lj6inv);
  kk_cut_lj3inv       = static_cast<KK_FLOAT>(cut_lj3inv);
  kk_cut_lj_inner6inv = static_cast<KK_FLOAT>(cut_lj_inner6inv);
  kk_cut_lj_inner3inv = static_cast<KK_FLOAT>(cut_lj_inner3inv);
  kk_g_ewald          = static_cast<KK_FLOAT>(g_ewald);

  Kokkos::deep_copy(d_cut_ljsq,   cut_ljsq);
  Kokkos::deep_copy(d_cut_coulsq, cut_coulsq);

  // error if rRESPA with inner levels

  if (update->whichflag == 1 &&
      utils::strmatch(update->integrate_style, "^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner  >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR, "Cannot use Kokkos pair style with rRESPA inner/middle");
  }

  // adjust neighbor list request for KOKKOS

  neighflag  = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCharmmfswCoulLongKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCharmmfswCoulLong::init_one(i, j);

  k_params.view_host()(i,j).lj1      = lj1[i][j];
  k_params.view_host()(i,j).lj2      = lj2[i][j];
  k_params.view_host()(i,j).lj3      = lj3[i][j];
  k_params.view_host()(i,j).lj4      = lj4[i][j];
  k_params.view_host()(i,j).cut_ljsq  = cut_ljsq;
  k_params.view_host()(i,j).cut_coulsq = cut_coulsq;

  k_params.view_host()(j,i) = k_params.view_host()(i,j);

  if (i < MAX_TYPES_STACKPARAMS+1 && j < MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j]      = m_params[j][i] = k_params.view_host()(i,j);
    m_cutsq[j][i]       = m_cutsq[i][j]  = cutone * cutone;
    m_cut_ljsq[j][i]    = m_cut_ljsq[i][j]    = cut_ljsq;
    m_cut_coulsq[j][i]  = m_cut_coulsq[i][j]  = cut_coulsq;
  }

  k_cutsq.view_host()(i,j) = k_cutsq.view_host()(j,i) = cutone * cutone;
  k_cutsq.modify_host();
  k_params.modify_host();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCharmmfswCoulLongKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCharmmfswCoulLongKokkos<LMPHostType>;
#endif
}
