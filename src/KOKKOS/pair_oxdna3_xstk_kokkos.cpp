/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_oxdna3_xstk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_npair_kokkos.h"
#include "fix_oxdna_prime_neighs_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;
using MathConst::MY_PI;

// NOTE: Regarding many math ops. I've added in alot of fma (Kokkos:fma)
// fused-multipy-add ops to reduce roundoffs and FP ops. More to the point,
// these were flagged up in Nsight Compute (with Source&SASS) as they really ate up
// Live Register Usage.
// fma (Kokkos::fma) fuses multiply-add operations: Kokkos::fma(x,y,z) = x*y + z,
// but with only one FP rounding error and one instruction instead of two.
// Additionally, we can use:
// sin^2(theta) = 1 - cos^2(theta) to rejig into fma format.
// Copilot (AI) was most definitly used to sweep through and apply the majoirty of
// these changes, vastly cutting down on my own errors/recompiling, as well as oversights
// I could easily miss - same applies to the bool terms I pull out for the sake of
// Live Register Usage.
// Local .sh regression tests have been ran (full FP64) to confirm numerics aren't altered.

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna3XstkKokkos<DeviceType>::PairOxdna3XstkKokkos(LAMMPS *lmp) : PairOxdna3Xstk(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  screened_pair_count = 0;
  fix_oxdna_lrfKK = nullptr;
  fix_oxdna_npairKK = nullptr;
  fix_oxdna_prime_neighsKK = nullptr;
  last_prime_neighs_xstk3_lastcall = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna3XstkKokkos<DeviceType>::~PairOxdna3XstkKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3XstkKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.template view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);

  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK | TORQUE_MASK);

  x = atomKK->k_x.template view<DeviceType>();
  f = atomKK->k_f.template view<DeviceType>();
  torque = atomKK->k_torque.template view<DeviceType>();
  type = atomKK->k_type.template view<DeviceType>();

  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  int need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    dup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterDuplicated>(f);
    dup_torque = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterDuplicated>(torque);
  } else {
    ndup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_torque = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(torque);
  }

  copymode = 1;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->k_nx.template view<DeviceType>();
  d_ny_xtrct = fix_oxdna_lrfKK->k_ny.template view<DeviceType>();
  d_nz_xtrct = fix_oxdna_lrfKK->k_nz.template view<DeviceType>();

  // Use the oxdna npair screened list on all backends.
  screened_pair_count = fix_oxdna_npairKK->screened_pair_count;
  d_pairs_screened = fix_oxdna_npairKK->k_pairs_screened.template view<DeviceType>();

  // Then get the precomputed 3'/5' neighbor map lookups for the screened npair list.
  // Done here (not in pre_force) so the pair's own list is always used,
  // ensuring ib-index correspondence between precompute and kernel.
  if (last_prime_neighs_xstk3_lastcall != neighbor->lastcall) {
    fix_oxdna_prime_neighsKK->compute_prime_neighs_oxdna3_xstk(list);
    last_prime_neighs_xstk3_lastcall = neighbor->lastcall;
    d_prime_neighs_oxdna3_xstk = fix_oxdna_prime_neighsKK->d_prime_neighs_oxdna3_xstk;
  }

  // loop over neighbors of my atoms for compute functors

  EV_FLOAT ev;

  // Launch from screened npair pairs regardless of backend.
  auto run_compute_screened = [&](auto screened_tag, auto evflag_tag) {
    constexpr int EVFLAG = decltype(evflag_tag)::value;
    if constexpr (EVFLAG) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, decltype(screened_tag)>(0,screened_pair_count),*this,ev);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, decltype(screened_tag)>(0,screened_pair_count),*this);
    }
  };

  auto run_compute_by_flags = [&](auto neighflag_tag, auto newtonpair_tag, auto evflag_tag) {
    constexpr int NEIGHFLAG = decltype(neighflag_tag)::value;
    constexpr int NEWTON_PAIR = decltype(newtonpair_tag)::value;
    constexpr int EVFLAG = decltype(evflag_tag)::value;
    run_compute_screened(TagPairOxdna3XstkComputeNpair<NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
  };

  const int dispatch_neigh =
      (neighflag == HALF) ? 0 :
      (neighflag == HALFTHREAD) ? 1 :
      (neighflag == FULL) ? 2 : -1;

  if (dispatch_neigh < 0) {
    error->all(FLERR, "Unsupported neighbor flag in pair oxdna3/xstk/kk");
  }

  const int dispatch_key = (evflag ? 8 : 0) | (newton_pair ? 4 : 0) | dispatch_neigh;
  switch (dispatch_key) {
    case 0: run_compute_by_flags(std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 1: run_compute_by_flags(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 2: run_compute_by_flags(std::integral_constant<int,FULL>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 4: run_compute_by_flags(std::integral_constant<int,HALF>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 5: run_compute_by_flags(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 6: run_compute_by_flags(std::integral_constant<int,FULL>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 8: run_compute_by_flags(std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 9: run_compute_by_flags(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 10: run_compute_by_flags(std::integral_constant<int,FULL>{},      std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 12: run_compute_by_flags(std::integral_constant<int,HALF>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    case 13: run_compute_by_flags(std::integral_constant<int,HALFTHREAD>{},std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    case 14: run_compute_by_flags(std::integral_constant<int,FULL>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    default: error->all(FLERR, "Internal dispatch error in pair oxdna3/xstk/kk");
  }

  if (need_dup) {
    Kokkos::Experimental::contribute(f, dup_f);
    Kokkos::Experimental::contribute(torque, dup_torque);
  }

  if (eflag_global) eng_vdwl += ev.evdwl;
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;

  // free duplicated memory
  if (need_dup) {
    dup_f        = decltype(dup_f)();
    dup_torque   = decltype(dup_torque)();
    dup_eatom    = decltype(dup_eatom)();
    dup_vatom    = decltype(dup_vatom)();
  }
}

/* ----------------------------------------------------------------------
   ComputeNpair Functor(s) - first we have all the xstk_* terms I don't
   really like....
   But these seem the best option performance-wise so far.
-------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_preradial_terms(
  KK_FLOAT &r_bsbs, KK_FLOAT &rinv_bsbs,
  KK_FLOAT (&delr_bsbs)[3], KK_FLOAT (&delr_bsbs_norm)[3],
  KK_FLOAT (&ra_cbs)[3], KK_FLOAT (&rb_cbs)[3],
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  const int &a, const int &b, const int &atype, const int &btype) const
{

  const KK_FLOAT dx_cbs_pur_oxdna3 = static_cast<KK_FLOAT>(0.43);
  const KK_FLOAT dx_cbs_pyr_oxdna3 = static_cast<KK_FLOAT>(0.37);

  const int anuc = atype % 4;
  const KK_FLOAT a_shift = (anuc == 0 || anuc == 2) ? dx_cbs_pyr_oxdna3 : dx_cbs_pur_oxdna3;
  ra_cbs[0] = a_shift * a_nx[0];
  ra_cbs[1] = a_shift * a_nx[1];
  ra_cbs[2] = a_shift * a_nx[2];

  const int bnuc = btype % 4;
  const KK_FLOAT b_shift = (bnuc == 0 || bnuc == 2) ? dx_cbs_pyr_oxdna3 : dx_cbs_pur_oxdna3;
  rb_cbs[0] = b_shift * b_nx[0];
  rb_cbs[1] = b_shift * b_nx[1];
  rb_cbs[2] = b_shift * b_nx[2];

  delr_bsbs[0] = x(a,0) + ra_cbs[0] - x(b,0) - rb_cbs[0];
  delr_bsbs[1] = x(a,1) + ra_cbs[1] - x(b,1) - rb_cbs[1];
  delr_bsbs[2] = x(a,2) + ra_cbs[2] - x(b,2) - rb_cbs[2];

  const KK_FLOAT rsq_bsbs = Kokkos::fma(delr_bsbs[2], delr_bsbs[2],
    Kokkos::fma(delr_bsbs[1], delr_bsbs[1], delr_bsbs[0] * delr_bsbs[0]));
  if (rsq_bsbs <= static_cast<KK_FLOAT>(0.0)) return false;

  rinv_bsbs = static_cast<KK_FLOAT>(1.0) / sqrtf(rsq_bsbs);
  r_bsbs = rsq_bsbs * rinv_bsbs;
  delr_bsbs_norm[0] = delr_bsbs[0] * rinv_bsbs;
  delr_bsbs_norm[1] = delr_bsbs[1] * rinv_bsbs;
  delr_bsbs_norm[2] = delr_bsbs[2] * rinv_bsbs;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_radial_terms(const int &atype, const int &btype,
  const int &a3ptype, const int &a5ptype,
  const int &b3ptype, const int &b5ptype,
  const KK_FLOAT &r_hb,
  KK_FLOAT &f2_33, KK_FLOAT &f2_55,
  KK_FLOAT &df2_33, KK_FLOAT &df2_55) const
{
  const auto& p_xstk = d_params_xstk(atype, btype);
  const auto& p_33 = d_params_33(a3ptype, atype, btype, b3ptype);
  const auto& p_55 = d_params_55(a5ptype, atype, btype, b5ptype);

  const KK_FLOAT l_k_xst = p_xstk.k_xst;
  const KK_FLOAT l_cut_xst_0 = p_33.cut_xst_0;
  const KK_FLOAT l_cut_xst_lc = p_33.cut_xst_lc;
  const KK_FLOAT l_cut_xst_hc = p_33.cut_xst_hc;
  const KK_FLOAT l_cut_xst_lo = p_33.cut_xst_lo;
  const KK_FLOAT l_cut_xst_hi = p_33.cut_xst_hi;
  const KK_FLOAT l_b_xst_lo = p_xstk.b_xst_lo;
  const KK_FLOAT l_b_xst_hi = p_xstk.b_xst_hi;
  const KK_FLOAT l_cut_xst_c = p_33.cut_xst_c;

  f2_33 = F2_KK(r_hb, l_k_xst, l_cut_xst_0, l_cut_xst_lc, l_cut_xst_hc,
                l_cut_xst_lo, l_cut_xst_hi, l_b_xst_lo, l_b_xst_hi, l_cut_xst_c, df2_33);

  f2_55 = F2_KK(r_hb, l_k_xst, p_55.cut_xst_0, p_55.cut_xst_lc, p_55.cut_xst_hc,
                p_55.cut_xst_lo, p_55.cut_xst_hi, l_b_xst_lo, l_b_xst_hi, p_55.cut_xst_c, df2_55);

  if (f2_33 == static_cast<KK_FLOAT>(0.0) && f2_55 == static_cast<KK_FLOAT>(0.0)) return false;

  // df2_33 = DF2_KK(r_hb, l_k_xst, l_cut_xst_0, l_cut_xst_lc, l_cut_xst_hc,
  //                 l_cut_xst_lo, l_cut_xst_hi, l_b_xst_lo, l_b_xst_hi);

  // df2_55 = DF2_KK(r_hb, l_k_xst, l_cut_xst_0, l_cut_xst_lc, p_55.cut_xst_hc,
  //                 p_55.cut_xst_lo, p_55.cut_xst_hi, l_b_xst_lo, l_b_xst_hi);

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta1_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  KK_FLOAT &f4t1, KK_FLOAT &df4t1) const
{
  KK_FLOAT cost1 = -Kokkos::fma(a_nx[2], b_nx[2], Kokkos::fma(a_nx[1], b_nx[1], a_nx[0] * b_nx[0]));
  if (cost1 > static_cast<KK_FLOAT>(1.0)) cost1 = static_cast<KK_FLOAT>(1.0);
  if (cost1 < static_cast<KK_FLOAT>(-1.0)) cost1 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta1 = acos(cost1);

  const auto& p_xstk = d_params_xstk(atype, btype);

  const KK_FLOAT l_a_xst1 = p_xstk.a_xst1;
  const KK_FLOAT l_theta_xst1_0 = p_xstk.theta_xst1_0;
  const KK_FLOAT l_dtheta_xst1_ast = p_xstk.dtheta_xst1_ast;
  const KK_FLOAT l_b_xst1 = p_xstk.b_xst1;
  const KK_FLOAT l_dtheta_xst1_c = p_xstk.dtheta_xst1_c;

  f4t1 = F4_KK(theta1, l_a_xst1, l_theta_xst1_0, l_dtheta_xst1_ast, l_b_xst1, l_dtheta_xst1_c, df4t1);
  if (f4t1 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin1_sq = Kokkos::fma(-cost1, cost1, static_cast<KK_FLOAT>(1.0));
  if (sin1_sq < static_cast<KK_FLOAT>(0.0)) sin1_sq = static_cast<KK_FLOAT>(0.0);
  if (sin1_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin1 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin1_sq));
  // df4t1 = DF4_KK(theta1, l_a_xst1, l_theta_xst1_0, l_dtheta_xst1_ast,
  //                l_b_xst1, l_dtheta_xst1_c) * rsin1;
  df4t1 *= rsin1;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta2_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &cost2, KK_FLOAT &f4t2, KK_FLOAT &df4t2) const
{
  cost2 = -Kokkos::fma(a_nx[2], delr_hb_norm[2], Kokkos::fma(a_nx[1], delr_hb_norm[1], a_nx[0] * delr_hb_norm[0]));
  if (cost2 > static_cast<KK_FLOAT>(1.0)) cost2 = static_cast<KK_FLOAT>(1.0);
  if (cost2 < static_cast<KK_FLOAT>(-1.0)) cost2 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta2 = acos(cost2);

  const auto& p_xstk = d_params_xstk(atype, btype);

  const KK_FLOAT l_a_xst2 = p_xstk.a_xst2;
  const KK_FLOAT l_theta_xst2_0 = p_xstk.theta_xst2_0;
  const KK_FLOAT l_dtheta_xst2_ast = p_xstk.dtheta_xst2_ast;
  const KK_FLOAT l_b_xst2 = p_xstk.b_xst2;
  const KK_FLOAT l_dtheta_xst2_c = p_xstk.dtheta_xst2_c;

  f4t2 = F4_KK(theta2, l_a_xst2, l_theta_xst2_0, l_dtheta_xst2_ast, l_b_xst2, l_dtheta_xst2_c, df4t2);
  if (f4t2 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin2_sq = Kokkos::fma(-cost2, cost2, static_cast<KK_FLOAT>(1.0));
  if (sin2_sq < static_cast<KK_FLOAT>(0.0)) sin2_sq = static_cast<KK_FLOAT>(0.0);
  if (sin2_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin2 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin2_sq));
  // df4t2 = DF4_KK(theta2, l_a_xst2, l_theta_xst2_0,
  //                l_dtheta_xst2_ast, l_b_xst2, l_dtheta_xst2_c) * rsin2;
  df4t2 *= rsin2;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta3_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &cost3, KK_FLOAT &f4t3, KK_FLOAT &df4t3) const
{
  cost3 = Kokkos::fma(b_nx[2], delr_hb_norm[2], Kokkos::fma(b_nx[1], delr_hb_norm[1], b_nx[0] * delr_hb_norm[0]));
  if (cost3 > static_cast<KK_FLOAT>(1.0)) cost3 = static_cast<KK_FLOAT>(1.0);
  if (cost3 < static_cast<KK_FLOAT>(-1.0)) cost3 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta3 = acos(cost3);

  const auto& p_xstk = d_params_xstk(atype, btype);

  const KK_FLOAT l_a_xst3 = p_xstk.a_xst3;
  const KK_FLOAT l_theta_xst3_0 = p_xstk.theta_xst3_0;
  const KK_FLOAT l_dtheta_xst3_ast = p_xstk.dtheta_xst3_ast;
  const KK_FLOAT l_b_xst3 = p_xstk.b_xst3;
  const KK_FLOAT l_dtheta_xst3_c = p_xstk.dtheta_xst3_c;

  f4t3 = F4_KK(theta3, l_a_xst3, l_theta_xst3_0, l_dtheta_xst3_ast, l_b_xst3, l_dtheta_xst3_c, df4t3);
  if (f4t3 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin3_sq = Kokkos::fma(-cost3, cost3, static_cast<KK_FLOAT>(1.0));
  if (sin3_sq < static_cast<KK_FLOAT>(0.0)) sin3_sq = static_cast<KK_FLOAT>(0.0);
  if (sin3_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin3 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin3_sq));
  // df4t3 = DF4_KK(theta3, l_a_xst3, l_theta_xst3_0,
  //                l_dtheta_xst3_ast, l_b_xst3, l_dtheta_xst3_c) * rsin3;
  df4t3 *= rsin3;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta4_terms(const int &atype, const int &btype,
  const int &a3ptype, const int &a5ptype,
  const int &b3ptype, const int &b5ptype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  KK_FLOAT &f4t4_33, KK_FLOAT &f4t4_55,
  KK_FLOAT &df4t4_33, KK_FLOAT &df4t4_55) const
{
  KK_FLOAT cost4 = Kokkos::fma(a_nz[2], b_nz[2], Kokkos::fma(a_nz[1], b_nz[1], a_nz[0] * b_nz[0]));
  if (cost4 > static_cast<KK_FLOAT>(1.0)) cost4 = static_cast<KK_FLOAT>(1.0);
  if (cost4 < static_cast<KK_FLOAT>(-1.0)) cost4 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta4 = acos(cost4);

  const auto& p_33 = d_params_33(a3ptype, atype, btype, b3ptype);
  const auto& p_55 = d_params_55(a5ptype, atype, btype, b5ptype);

  const KK_FLOAT l_a_xst4_33 = p_33.a_xst4;
  const KK_FLOAT l_theta_xst4_0_33 = p_33.theta_xst4_0;
  const KK_FLOAT l_dtheta_xst4_ast_33 = p_33.dtheta_xst4_ast;
  const KK_FLOAT l_b_xst4_33 = p_33.b_xst4;
  const KK_FLOAT l_dtheta_xst4_c_33 = p_33.dtheta_xst4_c;
  f4t4_33 = F4_KK(theta4, l_a_xst4_33, l_theta_xst4_0_33, l_dtheta_xst4_ast_33, l_b_xst4_33, l_dtheta_xst4_c_33, df4t4_33);

  const KK_FLOAT l_a_xst4_55 = p_55.a_xst4;
  const KK_FLOAT l_theta_xst4_0_55 = p_55.theta_xst4_0;
  const KK_FLOAT l_dtheta_xst4_ast_55 = p_55.dtheta_xst4_ast;
  const KK_FLOAT l_b_xst4_55 = p_55.b_xst4;
  const KK_FLOAT l_dtheta_xst4_c_55 = p_55.dtheta_xst4_c;
  f4t4_55 = F4_KK(theta4, l_a_xst4_55, l_theta_xst4_0_55, l_dtheta_xst4_ast_55, l_b_xst4_55, l_dtheta_xst4_c_55, df4t4_55);
  if (f4t4_33 == static_cast<KK_FLOAT>(0.0) && f4t4_55 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin4_sq = Kokkos::fma(-cost4, cost4, static_cast<KK_FLOAT>(1.0));
  if (sin4_sq < static_cast<KK_FLOAT>(0.0)) sin4_sq = static_cast<KK_FLOAT>(0.0);
  if (sin4_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin4 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin4_sq));

  // df4t4_33 = DF4_KK(theta4, l_a_xst4_33, l_theta_xst4_0_33, l_dtheta_xst4_ast_33, l_b_xst4_33, l_dtheta_xst4_c_33) * rsin4;
  // df4t4_55 = DF4_KK(theta4, l_a_xst4_55, l_theta_xst4_0_55, l_dtheta_xst4_ast_55, l_b_xst4_55, l_dtheta_xst4_c_55) * rsin4;
  df4t4_33 *= rsin4;
  df4t4_55 *= rsin4;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta7_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &cost7,
  KK_FLOAT &f4t7_33, KK_FLOAT &f4t7_55,
  KK_FLOAT &df4t7_33, KK_FLOAT &df4t7_55) const
{
  cost7 = -Kokkos::fma(a_nz[2], delr_hb_norm[2], Kokkos::fma(a_nz[1], delr_hb_norm[1], a_nz[0] * delr_hb_norm[0]));
  if (cost7 > static_cast<KK_FLOAT>(1.0)) cost7 = static_cast<KK_FLOAT>(1.0);
  if (cost7 < static_cast<KK_FLOAT>(-1.0)) cost7 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta7 = acos(cost7);

  const auto& p_t7 = d_params_t7(atype, btype);

  const KK_FLOAT l_a_xst7 = p_t7.a_xst7;
  const KK_FLOAT l_theta_xst7_0_33 = p_t7.theta_xst7_0_33;
  const KK_FLOAT l_dtheta_xst7_ast = p_t7.dtheta_xst7_ast;
  const KK_FLOAT l_b_xst7 = p_t7.b_xst7;
  const KK_FLOAT l_dtheta_xst7_c = p_t7.dtheta_xst7_c;
  f4t7_33 = F4_KK(theta7, l_a_xst7, l_theta_xst7_0_33, l_dtheta_xst7_ast, l_b_xst7, l_dtheta_xst7_c, df4t7_33);
  const KK_FLOAT l_theta_xst7_0_55 = p_t7.theta_xst7_0_55;
  f4t7_55 = F4_KK(theta7, l_a_xst7, l_theta_xst7_0_55, l_dtheta_xst7_ast, l_b_xst7, l_dtheta_xst7_c, df4t7_55);
  if (f4t7_33 == static_cast<KK_FLOAT>(0.0) && f4t7_55 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin7_sq = Kokkos::fma(-cost7, cost7, static_cast<KK_FLOAT>(1.0));
  if (sin7_sq < static_cast<KK_FLOAT>(0.0)) sin7_sq = static_cast<KK_FLOAT>(0.0);
  if (sin7_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin7 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin7_sq));

  // df4t7_33 = DF4_KK(theta7, l_a_xst7, l_theta_xst7_0_33, l_dtheta_xst7_ast, l_b_xst7, l_dtheta_xst7_c) * rsin7;
  // df4t7_55 = DF4_KK(theta7, l_a_xst7, l_theta_xst7_0_55, l_dtheta_xst7_ast, l_b_xst7, l_dtheta_xst7_c) * rsin7;
  df4t7_33 *= rsin7;
  df4t7_55 *= rsin7;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdna3XstkKokkos<DeviceType>::xstk_theta8_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &cost8,
  KK_FLOAT &f4t8_33, KK_FLOAT &f4t8_55,
  KK_FLOAT &df4t8_33, KK_FLOAT &df4t8_55) const
{
  cost8 = Kokkos::fma(b_nz[2], delr_hb_norm[2], Kokkos::fma(b_nz[1], delr_hb_norm[1], b_nz[0] * delr_hb_norm[0]));
  if (cost8 > static_cast<KK_FLOAT>(1.0)) cost8 = static_cast<KK_FLOAT>(1.0);
  if (cost8 < static_cast<KK_FLOAT>(-1.0)) cost8 = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT theta8 = acos(cost8);

  const auto& p_t8 = d_params_t8(atype, btype);

  const KK_FLOAT l_a_xst8 = p_t8.a_xst8;
  const KK_FLOAT l_theta_xst8_0_33 = p_t8.theta_xst8_0_33;
  const KK_FLOAT l_dtheta_xst8_ast = p_t8.dtheta_xst8_ast;
  const KK_FLOAT l_b_xst8 = p_t8.b_xst8;
  const KK_FLOAT l_dtheta_xst8_c = p_t8.dtheta_xst8_c;
  f4t8_33 = F4_KK(theta8, l_a_xst8, l_theta_xst8_0_33, l_dtheta_xst8_ast, l_b_xst8, l_dtheta_xst8_c, df4t8_33);
  const KK_FLOAT l_theta_xst8_0_55 = p_t8.theta_xst8_0_55;
  f4t8_55 = F4_KK(theta8, l_a_xst8, l_theta_xst8_0_55, l_dtheta_xst8_ast, l_b_xst8, l_dtheta_xst8_c, df4t8_55);
  if (f4t8_33 == static_cast<KK_FLOAT>(0.0) && f4t8_55 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin8_sq = Kokkos::fma(-cost8, cost8, static_cast<KK_FLOAT>(1.0));
  if (sin8_sq < static_cast<KK_FLOAT>(0.0)) sin8_sq = static_cast<KK_FLOAT>(0.0);
  if (sin8_sq <= static_cast<KK_FLOAT>(1.0e-12)) return false;
  const KK_FLOAT rsin8 = static_cast<KK_FLOAT>(Kokkos::rsqrt(sin8_sq));

  // df4t8_33 = DF4_KK(theta8, l_a_xst8, l_theta_xst8_0_33, l_dtheta_xst8_ast, l_b_xst8, l_dtheta_xst8_c) * rsin8;
  // df4t8_55 = DF4_KK(theta8, l_a_xst8, l_theta_xst8_0_55, l_dtheta_xst8_ast, l_b_xst8, l_dtheta_xst8_c) * rsin8;
  df4t8_33 *= rsin8;
  df4t8_55 *= rsin8;

  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdna3XstkKokkos<DeviceType>::xstk_force_contrib(const KK_FLOAT &f2_33, const KK_FLOAT &f2_55,
  const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
  const KK_FLOAT &f4t4_33, const KK_FLOAT &f4t4_55,
  const KK_FLOAT &f4t7_33, const KK_FLOAT &f4t7_55,
  const KK_FLOAT &f4t8_33, const KK_FLOAT &f4t8_55,
  const KK_FLOAT &df2_33, const KK_FLOAT &df2_55,
  const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
  const KK_FLOAT &df4t7_33, const KK_FLOAT &df4t7_55,
  const KK_FLOAT &df4t8_33, const KK_FLOAT &df4t8_55,
  const KK_FLOAT &rinv_hb, const KK_FLOAT &factor_lj,
  const KK_FLOAT &cost2, const KK_FLOAT &cost3, const KK_FLOAT &cost7, const KK_FLOAT &cost8,
  const KK_FLOAT (&delr_hb)[3], const KK_FLOAT (&delr_hb_norm)[3],
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&ra_chb)[3], const KK_FLOAT (&rb_chb)[3],
  KK_ACC_FLOAT (&delf)[3], KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  const KK_FLOAT mixsum = f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55;
  const KK_FLOAT rsum_df2 = df2_33 * f4t4_33 * f4t7_33 * f4t8_33 + df2_55 * f4t4_55 * f4t7_55 * f4t8_55;

  // NOTE: Our early rejection in previous functions already covers the "if !theta" checks here that
  // are seen in vanilla.
  KK_FLOAT finc = -f4t1 * f4t2 * f4t3 * rsum_df2 * rinv_hb * factor_lj;
  delf[0] = Kokkos::fma(delr_hb[0], finc, delf[0]);
  delf[1] = Kokkos::fma(delr_hb[1], finc, delf[1]);
  delf[2] = Kokkos::fma(delr_hb[2], finc, delf[2]);

  finc = -f4t1 * df4t2 * f4t3 * mixsum * rinv_hb * factor_lj;
  delf[0] = Kokkos::fma(Kokkos::fma(delr_hb_norm[0], cost2, a_nx[0]), finc, delf[0]);
  delf[1] = Kokkos::fma(Kokkos::fma(delr_hb_norm[1], cost2, a_nx[1]), finc, delf[1]);
  delf[2] = Kokkos::fma(Kokkos::fma(delr_hb_norm[2], cost2, a_nx[2]), finc, delf[2]);

  finc = -f4t1 * f4t2 * df4t3 * mixsum * rinv_hb * factor_lj;
  delf[0] = Kokkos::fma(Kokkos::fma(delr_hb_norm[0], cost3, -b_nx[0]), finc, delf[0]);
  delf[1] = Kokkos::fma(Kokkos::fma(delr_hb_norm[1], cost3, -b_nx[1]), finc, delf[1]);
  delf[2] = Kokkos::fma(Kokkos::fma(delr_hb_norm[2], cost3, -b_nx[2]), finc, delf[2]);

  const KK_FLOAT t7sum = f2_33 * f4t4_33 * df4t7_33 * f4t8_33 + f2_55 * f4t4_55 * df4t7_55 * f4t8_55;
  finc = -f4t1 * f4t2 * f4t3 * t7sum * rinv_hb * factor_lj;
  delf[0] = Kokkos::fma(Kokkos::fma(delr_hb_norm[0], cost7, a_nz[0]), finc, delf[0]);
  delf[1] = Kokkos::fma(Kokkos::fma(delr_hb_norm[1], cost7, a_nz[1]), finc, delf[1]);
  delf[2] = Kokkos::fma(Kokkos::fma(delr_hb_norm[2], cost7, a_nz[2]), finc, delf[2]);

  const KK_FLOAT t8sum = f2_33 * f4t4_33 * f4t7_33 * df4t8_33 + f2_55 * f4t4_55 * f4t7_55 * df4t8_55;
  finc = -f4t1 * f4t2 * f4t3 * t8sum * rinv_hb * factor_lj;
  delf[0] = Kokkos::fma(Kokkos::fma(delr_hb_norm[0], cost8, -b_nz[0]), finc, delf[0]);
  delf[1] = Kokkos::fma(Kokkos::fma(delr_hb_norm[1], cost8, -b_nz[1]), finc, delf[1]);
  delf[2] = Kokkos::fma(Kokkos::fma(delr_hb_norm[2], cost8, -b_nz[2]), finc, delf[2]);

  delta[0] = Kokkos::fma(ra_chb[1], delf[2], -ra_chb[2] * delf[1]);
  delta[1] = Kokkos::fma(ra_chb[2], delf[0], -ra_chb[0] * delf[2]);
  delta[2] = Kokkos::fma(ra_chb[0], delf[1], -ra_chb[1] * delf[0]);

  deltb[0] = Kokkos::fma(rb_chb[1], delf[2], -rb_chb[2] * delf[1]);
  deltb[1] = Kokkos::fma(rb_chb[2], delf[0], -rb_chb[0] * delf[2]);
  deltb[2] = Kokkos::fma(rb_chb[0], delf[1], -rb_chb[1] * delf[0]);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdna3XstkKokkos<DeviceType>::xstk_torque_contrib(const KK_FLOAT &f2_33, const KK_FLOAT &f2_55,
  const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
  const KK_FLOAT &f4t4_33, const KK_FLOAT &f4t4_55,
  const KK_FLOAT &f4t7_33, const KK_FLOAT &f4t7_55,
  const KK_FLOAT &f4t8_33, const KK_FLOAT &f4t8_55,
  const KK_FLOAT &df4t1, const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
  const KK_FLOAT &df4t4_33, const KK_FLOAT &df4t4_55,
  const KK_FLOAT &df4t7_33, const KK_FLOAT &df4t7_55,
  const KK_FLOAT &df4t8_33, const KK_FLOAT &df4t8_55,
  const KK_FLOAT &factor_lj,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&delr_hb_norm)[3],
  KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  delta[0] = static_cast<KK_ACC_FLOAT>(0.0);
  delta[1] = static_cast<KK_ACC_FLOAT>(0.0);
  delta[2] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[0] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[1] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[2] = static_cast<KK_ACC_FLOAT>(0.0);

  const KK_FLOAT mixsum = f2_33 * f4t4_33 * f4t7_33 * f4t8_33 + f2_55 * f4t4_55 * f4t7_55 * f4t8_55;
  const KK_FLOAT pure_torque_common = f4t1 * f4t2 * f4t3 * factor_lj;

  // NOTE: Our early rejection in previous functions already covers the "if !theta" checks here that
  // are seen in vanilla.
  const KK_FLOAT tpair1 = -df4t1 * f4t2 * f4t3 * mixsum * factor_lj;
  const KK_FLOAT t1dir0 = Kokkos::fma(a_nx[1], b_nx[2], -a_nx[2] * b_nx[1]);
  const KK_FLOAT t1dir1 = Kokkos::fma(a_nx[2], b_nx[0], -a_nx[0] * b_nx[2]);
  const KK_FLOAT t1dir2 = Kokkos::fma(a_nx[0], b_nx[1], -a_nx[1] * b_nx[0]);
  delta[0] += t1dir0 * tpair1;
  delta[1] += t1dir1 * tpair1;
  delta[2] += t1dir2 * tpair1;
  deltb[0] += t1dir0 * tpair1;
  deltb[1] += t1dir1 * tpair1;
  deltb[2] += t1dir2 * tpair1;

  const KK_FLOAT tpair2 = -f4t1 * df4t2 * f4t3 * mixsum * factor_lj;
  const KK_FLOAT t2dir0 = Kokkos::fma(a_nx[1], delr_hb_norm[2], -a_nx[2] * delr_hb_norm[1]);
  const KK_FLOAT t2dir1 = Kokkos::fma(a_nx[2], delr_hb_norm[0], -a_nx[0] * delr_hb_norm[2]);
  const KK_FLOAT t2dir2 = Kokkos::fma(a_nx[0], delr_hb_norm[1], -a_nx[1] * delr_hb_norm[0]);
  delta[0] += t2dir0 * tpair2;
  delta[1] += t2dir1 * tpair2;
  delta[2] += t2dir2 * tpair2;

  const KK_FLOAT tpair3 = -f4t1 * f4t2 * df4t3 * mixsum * factor_lj;
  const KK_FLOAT t3dir0 = Kokkos::fma(b_nx[1], delr_hb_norm[2], -b_nx[2] * delr_hb_norm[1]);
  const KK_FLOAT t3dir1 = Kokkos::fma(b_nx[2], delr_hb_norm[0], -b_nx[0] * delr_hb_norm[2]);
  const KK_FLOAT t3dir2 = Kokkos::fma(b_nx[0], delr_hb_norm[1], -b_nx[1] * delr_hb_norm[0]);
  deltb[0] += t3dir0 * tpair3;
  deltb[1] += t3dir1 * tpair3;
  deltb[2] += t3dir2 * tpair3;

  const KK_FLOAT t4sum = f2_33 * df4t4_33 * f4t7_33 * f4t8_33 + f2_55 * df4t4_55 * f4t7_55 * f4t8_55;
  const KK_FLOAT tpair4 = -pure_torque_common * t4sum;
  const KK_FLOAT t4dir0 = Kokkos::fma(b_nz[1], a_nz[2], -b_nz[2] * a_nz[1]);
  const KK_FLOAT t4dir1 = Kokkos::fma(b_nz[2], a_nz[0], -b_nz[0] * a_nz[2]);
  const KK_FLOAT t4dir2 = Kokkos::fma(b_nz[0], a_nz[1], -b_nz[1] * a_nz[0]);
  delta[0] += t4dir0 * tpair4;
  delta[1] += t4dir1 * tpair4;
  delta[2] += t4dir2 * tpair4;
  deltb[0] += t4dir0 * tpair4;
  deltb[1] += t4dir1 * tpair4;
  deltb[2] += t4dir2 * tpair4;

  const KK_FLOAT t7sum = f2_33 * f4t4_33 * df4t7_33 * f4t8_33 + f2_55 * f4t4_55 * df4t7_55 * f4t8_55;
  const KK_FLOAT tpair7 = -pure_torque_common * t7sum;
  const KK_FLOAT t7dir0 = Kokkos::fma(a_nz[1], delr_hb_norm[2], -a_nz[2] * delr_hb_norm[1]);
  const KK_FLOAT t7dir1 = Kokkos::fma(a_nz[2], delr_hb_norm[0], -a_nz[0] * delr_hb_norm[2]);
  const KK_FLOAT t7dir2 = Kokkos::fma(a_nz[0], delr_hb_norm[1], -a_nz[1] * delr_hb_norm[0]);
  delta[0] += t7dir0 * tpair7;
  delta[1] += t7dir1 * tpair7;
  delta[2] += t7dir2 * tpair7;

  const KK_FLOAT t8sum = f2_33 * f4t4_33 * f4t7_33 * df4t8_33 + f2_55 * f4t4_55 * f4t7_55 * df4t8_55;
  const KK_FLOAT tpair8 = -pure_torque_common * t8sum;
  const KK_FLOAT t8dir0 = Kokkos::fma(b_nz[1], delr_hb_norm[2], -b_nz[2] * delr_hb_norm[1]);
  const KK_FLOAT t8dir1 = Kokkos::fma(b_nz[2], delr_hb_norm[0], -b_nz[0] * delr_hb_norm[2]);
  const KK_FLOAT t8dir2 = Kokkos::fma(b_nz[0], delr_hb_norm[1], -b_nz[1] * delr_hb_norm[0]);
  deltb[0] += t8dir0 * tpair8;
  deltb[1] += t8dir1 * tpair8;
  deltb[2] += t8dir2 * tpair8;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna3XstkKokkos<DeviceType>::operator()(TagPairOxdna3XstkComputeNpair<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair, EV_FLOAT &ev) const
{
  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const uint64_t pair = d_pairs_screened(ipair);
  const int a = static_cast<int>(pair >> 32);
  const int atype = type(a);
  int braw = static_cast<int>(pair & 0xffffffffu);
  const KK_FLOAT factor_lj = static_cast<KK_FLOAT>(special_lj[sbmask(braw)]);
  if (factor_lj == static_cast<KK_FLOAT>(0.0)) return;
  const int b = braw & NEIGHMASK;
  const int btype = type(b);

  const int a3idx = d_prime_neighs_oxdna3_xstk(ipair,0);
  const int a5idx = d_prime_neighs_oxdna3_xstk(ipair,1);
  const int b3idx = d_prime_neighs_oxdna3_xstk(ipair,2);
  const int b5idx = d_prime_neighs_oxdna3_xstk(ipair,3);

  const int a3ptype = (a3idx >= 0) ? type(a3idx) : 0;
  const int a5ptype = (a5idx >= 0) ? type(a5idx) : 0;
  const int b3ptype = (b3idx >= 0) ? type(b3idx) : 0;
  const int b5ptype = (b5idx >= 0) ? type(b5idx) : 0;

  KK_FLOAT a_nx[3], b_nx[3];
  a_nx[0] = d_nx_xtrct(a,0);
  a_nx[1] = d_nx_xtrct(a,1);
  a_nx[2] = d_nx_xtrct(a,2);
  b_nx[0] = d_nx_xtrct(b,0);
  b_nx[1] = d_nx_xtrct(b,1);
  b_nx[2] = d_nx_xtrct(b,2);

  KK_FLOAT r_bsbs, rinv_bsbs;
  KK_FLOAT ra_cbs[3], rb_cbs[3], delr_bsbs[3],delr_bsbs_norm[3];
  if (!xstk_preradial_terms(r_bsbs, rinv_bsbs, delr_bsbs, delr_bsbs_norm, ra_cbs, rb_cbs,
      a_nx, b_nx, a, b, atype, btype)) return;

  KK_FLOAT f2_33, f2_55, df2_33, df2_55;
  if (!xstk_radial_terms(atype, btype, a3ptype, a5ptype, b3ptype, b5ptype,
      r_bsbs, f2_33, f2_55, df2_33, df2_55)) return;

  KK_FLOAT f4t1, df4t1;
  if (!xstk_theta1_terms(atype, btype, a_nx, b_nx, f4t1, df4t1)) return;

  KK_FLOAT cost2, f4t2, df4t2;
  if (!xstk_theta2_terms(atype, btype, a_nx, delr_bsbs_norm, cost2, f4t2, df4t2)) return;

  KK_FLOAT cost3, f4t3, df4t3;
  if (!xstk_theta3_terms(atype, btype, b_nx, delr_bsbs_norm, cost3, f4t3, df4t3)) return;

  KK_FLOAT a_nz[3], b_nz[3];
  a_nz[0] = d_nz_xtrct(a,0);
  a_nz[1] = d_nz_xtrct(a,1);
  a_nz[2] = d_nz_xtrct(a,2);
  b_nz[0] = d_nz_xtrct(b,0);
  b_nz[1] = d_nz_xtrct(b,1);
  b_nz[2] = d_nz_xtrct(b,2);

  KK_FLOAT f4t4_33, f4t4_55, df4t4_33, df4t4_55;
  if (!xstk_theta4_terms(atype, btype, a3ptype, a5ptype, b3ptype, b5ptype,
      a_nz, b_nz, f4t4_33, f4t4_55, df4t4_33, df4t4_55)) return;

  KK_FLOAT cost7, f4t7_33, f4t7_55, df4t7_33, df4t7_55;
  if (!xstk_theta7_terms(atype, btype, a_nz, delr_bsbs_norm,
      cost7, f4t7_33, f4t7_55, df4t7_33, df4t7_55)) return;

  KK_FLOAT cost8, f4t8_33, f4t8_55, df4t8_33, df4t8_55;
  if (!xstk_theta8_terms(atype, btype, b_nz, delr_bsbs_norm,
      cost8, f4t8_33, f4t8_55, df4t8_33, df4t8_55)) return;

  const KK_FLOAT sum33 = f2_33 * f4t4_33 * f4t7_33 * f4t8_33;
  const KK_FLOAT sum55 = f2_55 * f4t4_55 * f4t7_55 * f4t8_55;
  const KK_FLOAT mixsum = sum33 + sum55;
  const KK_FLOAT evdwl = f4t1 * f4t2 * f4t3 * mixsum * factor_lj;
  if (evdwl == static_cast<KK_FLOAT>(0.0)) return;

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];
  delf[0] = static_cast<KK_ACC_FLOAT>(0.0);
  delf[1] = static_cast<KK_ACC_FLOAT>(0.0);
  delf[2] = static_cast<KK_ACC_FLOAT>(0.0);
  delta[0] = static_cast<KK_ACC_FLOAT>(0.0);
  delta[1] = static_cast<KK_ACC_FLOAT>(0.0);
  delta[2] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[0] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[1] = static_cast<KK_ACC_FLOAT>(0.0);
  deltb[2] = static_cast<KK_ACC_FLOAT>(0.0);

  xstk_force_contrib(
    f2_33, f2_55,
    f4t1, f4t2, f4t3,
    f4t4_33, f4t4_55,
    f4t7_33, f4t7_55,
    f4t8_33, f4t8_55,
    df2_33, df2_55,
    df4t2, df4t3,
    df4t7_33, df4t7_55,
    df4t8_33, df4t8_55,
    rinv_bsbs, factor_lj,
    cost2, cost3, cost7, cost8,
    delr_bsbs, delr_bsbs_norm,
    a_nx, b_nx, a_nz, b_nz,
    ra_cbs, rb_cbs,
    delf, delta, deltb);

  a_f(a,0) += delf[0];
  a_f(a,1) += delf[1];
  a_f(a,2) += delf[2];
  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];

  const bool do_newton_b = (NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || b < nlocal);
  if (do_newton_b) {
    a_f(b,0) -= delf[0];
    a_f(b,1) -= delf[1];
    a_f(b,2) -= delf[2];
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }

  if (EVFLAG) {
    if (eflag) {
      ev.evdwl += (do_newton_b ? static_cast<KK_ACC_FLOAT>(1.0) : static_cast<KK_ACC_FLOAT>(0.5)) * evdwl;
    }
    if (vflag_either || eflag_atom) {
      this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl,
        delf[0],delf[1],delf[2],x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
    }
  }

  xstk_torque_contrib(
    f2_33, f2_55,
    f4t1, f4t2, f4t3,
    f4t4_33, f4t4_55,
    f4t7_33, f4t7_55,
    f4t8_33, f4t8_55,
    df4t1, df4t2, df4t3,
    df4t4_33, df4t4_55,
    df4t7_33, df4t7_55,
    df4t8_33, df4t8_55,
    factor_lj,
    a_nx, b_nx, a_nz, b_nz,
    delr_bsbs_norm,
    delta, deltb);

  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];
  if (do_newton_b) {
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna3XstkKokkos<DeviceType>::operator()(TagPairOxdna3XstkComputeNpair<NEIGHFLAG,NEWTON_PAIR,EVFLAG>,
  const int &ipair) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>
    (TagPairOxdna3XstkComputeNpair<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ipair,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3XstkKokkos<DeviceType>::allocate()
{
  PairOxdna3Xstk::allocate();

  int n = atom->ntypes;

  // Allocate consolidated parameter views
  k_params_xstk = decltype(k_params_xstk)("PairOxdna3Xstk:params_xstk", n + 1, n + 1);
  k_params_33 = decltype(k_params_33)("PairOxdna3Xstk:params_33", n + 1, n + 1, n + 1, n + 1);
  k_params_55 = decltype(k_params_55)("PairOxdna3Xstk:params_55", n + 1, n + 1, n + 1, n + 1);
  k_params_t7 = decltype(k_params_t7)("PairOxdna3Xstk:params_t7", n + 1, n + 1);
  k_params_t8 = decltype(k_params_t8)("PairOxdna3Xstk:params_t8", n + 1, n + 1);

  d_params_xstk = k_params_xstk.template view<DeviceType>();
  d_params_33 = k_params_33.template view<DeviceType>();
  d_params_55 = k_params_55.template view<DeviceType>();
  d_params_t7 = k_params_t7.template view<DeviceType>();
  d_params_t8 = k_params_t8.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3XstkKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3XstkKokkos<DeviceType>::init_style()
{
  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);

  fix_oxdna_npairKK = nullptr;
  auto npair_fixes = modify->get_fix_by_style("^OXDNA/NPAIR/kk");
  if (npair_fixes.size() == 0) {
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(modify->add_fix("npair_kk all OXDNA/NPAIR/kk"));
  } else {
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(npair_fixes[0]);
  }
  if (!fix_oxdna_npairKK) error->all(FLERR, "Fix OXDNA/NPAIR/kk lookup failed");

  auto prime_fixes = modify->get_fix_by_style("^OXDNA/PRIME_NEIGHS/kk");
  if (prime_fixes.size() == 0) {
    fix_oxdna_prime_neighsKK =
      dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(modify->add_fix("prime_neighs_kk all OXDNA/PRIME_NEIGHS/kk"));
  } else {
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(prime_fixes[0]);
  }
  if (!fix_oxdna_prime_neighsKK) error->all(FLERR, "Fix OXDNA/PRIME_NEIGHS/kk lookup failed");

  // oxdna3/xstk always uses the npair screened list; force rebuilds on all backends.
  fix_oxdna_npairKK->set_force_screening_all_backends(true);
}

/* ----------------------------------------------------------------------
   All non-tetramer Kokkos views are set here within ::init_one, and
   the tetramer Kokkos views are set within ::coeff
------------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdna3XstkKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdna3Xstk::init_one(i,j);

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  // Populate the 2D parameter struct
  auto h_params_xstk = k_params_xstk.view_host();
  h_params_xstk(i, j).k_xst = k_xst[i][j]; h_params_xstk(j, i).k_xst = k_xst[j][i];
  h_params_xstk(i, j).b_xst_lo = b_xst_lo[i][j]; h_params_xstk(j, i).b_xst_lo = b_xst_lo[j][i];
  h_params_xstk(i, j).b_xst_hi = b_xst_hi[i][j]; h_params_xstk(j, i).b_xst_hi = b_xst_hi[j][i];

  h_params_xstk(i, j).a_xst1 = a_xst1[i][j]; h_params_xstk(j, i).a_xst1 = a_xst1[j][i];
  h_params_xstk(i, j).theta_xst1_0 = theta_xst1_0[i][j]; h_params_xstk(j, i).theta_xst1_0 = theta_xst1_0[j][i];
  h_params_xstk(i, j).dtheta_xst1_ast = dtheta_xst1_ast[i][j]; h_params_xstk(j, i).dtheta_xst1_ast = dtheta_xst1_ast[j][i];
  h_params_xstk(i, j).b_xst1 = b_xst1[i][j]; h_params_xstk(j, i).b_xst1 = b_xst1[j][i];
  h_params_xstk(i, j).dtheta_xst1_c = dtheta_xst1_c[i][j]; h_params_xstk(j, i).dtheta_xst1_c = dtheta_xst1_c[j][i];
  h_params_xstk(i, j).a_xst2 = a_xst2[i][j]; h_params_xstk(j, i).a_xst2 = a_xst2[j][i];
  h_params_xstk(i, j).theta_xst2_0 = theta_xst2_0[i][j]; h_params_xstk(j, i).theta_xst2_0 = theta_xst2_0[j][i];
  h_params_xstk(i, j).dtheta_xst2_ast = dtheta_xst2_ast[i][j]; h_params_xstk(j, i).dtheta_xst2_ast = dtheta_xst2_ast[j][i];
  h_params_xstk(i, j).b_xst2 = b_xst2[i][j]; h_params_xstk(j, i).b_xst2 = b_xst2[j][i];
  h_params_xstk(i, j).dtheta_xst2_c = dtheta_xst2_c[i][j]; h_params_xstk(j, i).dtheta_xst2_c = dtheta_xst2_c[j][i];
  h_params_xstk(i, j).a_xst3 = a_xst3[i][j]; h_params_xstk(j, i).a_xst3 = a_xst3[j][i];
  h_params_xstk(i, j).theta_xst3_0 = theta_xst3_0[i][j]; h_params_xstk(j, i).theta_xst3_0 = theta_xst3_0[j][i];
  h_params_xstk(i, j).dtheta_xst3_ast = dtheta_xst3_ast[i][j]; h_params_xstk(j, i).dtheta_xst3_ast = dtheta_xst3_ast[j][i];
  h_params_xstk(i, j).b_xst3 = b_xst3[i][j]; h_params_xstk(j, i).b_xst3 = b_xst3[j][i];
  h_params_xstk(i, j).dtheta_xst3_c = dtheta_xst3_c[i][j]; h_params_xstk(j, i).dtheta_xst3_c = dtheta_xst3_c[j][i];

  // Populate theta7 and theta8 parameter structs
  auto h_params_t7 = k_params_t7.view_host();
  h_params_t7(i, j).a_xst7 = a_xst7[i][j]; h_params_t7(j, i).a_xst7 = a_xst7[j][i];
  h_params_t7(i, j).theta_xst7_0_33 = theta_xst7_0_33[i][j]; h_params_t7(j, i).theta_xst7_0_33 = theta_xst7_0_33[j][i];
  h_params_t7(i, j).theta_xst7_0_55 = theta_xst7_0_55[i][j]; h_params_t7(j, i).theta_xst7_0_55 = theta_xst7_0_55[j][i];
  h_params_t7(i, j).dtheta_xst7_ast = dtheta_xst7_ast[i][j]; h_params_t7(j, i).dtheta_xst7_ast = dtheta_xst7_ast[j][i];
  h_params_t7(i, j).b_xst7 = b_xst7[i][j]; h_params_t7(j, i).b_xst7 = b_xst7[j][i];
  h_params_t7(i, j).dtheta_xst7_c = dtheta_xst7_c[i][j]; h_params_t7(j, i).dtheta_xst7_c = dtheta_xst7_c[j][i];

  auto h_params_t8 = k_params_t8.view_host();
  h_params_t8(i, j).a_xst8 = a_xst8[i][j]; h_params_t8(j, i).a_xst8 = a_xst8[j][i];
  h_params_t8(i, j).theta_xst8_0_33 = theta_xst8_0_33[i][j]; h_params_t8(j, i).theta_xst8_0_33 = theta_xst8_0_33[j][i];
  h_params_t8(i, j).theta_xst8_0_55 = theta_xst8_0_55[i][j]; h_params_t8(j, i).theta_xst8_0_55 = theta_xst8_0_55[j][i];
  h_params_t8(i, j).dtheta_xst8_ast = dtheta_xst8_ast[i][j]; h_params_t8(j, i).dtheta_xst8_ast = dtheta_xst8_ast[j][i];
  h_params_t8(i, j).b_xst8 = b_xst8[i][j]; h_params_t8(j, i).b_xst8 = b_xst8[j][i];
  h_params_t8(i, j).dtheta_xst8_c = dtheta_xst8_c[i][j]; h_params_t8(j, i).dtheta_xst8_c = dtheta_xst8_c[j][i];

  k_params_xstk.template modify<LMPHostType>();
  k_params_t7.template modify<LMPHostType>();
  k_params_t8.template modify<LMPHostType>();

  k_params_xstk.template sync<DeviceType>();
  k_params_t7.template sync<DeviceType>();
  k_params_t8.template sync<DeviceType>();

  return cutone;
}

/* ----------------------------------------------------------------------
   All tetramer Kokkos views are set here within ::coeff,
   and the non-tetramer Kokkos views are set within ::init_one
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna3XstkKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairOxdna3Xstk::coeff(narg,arg);

  const int n = atom->ntypes;

  // Populate the 4D parameter structs
  auto h_params_33 = k_params_33.view_host();
  auto h_params_55 = k_params_55.view_host();
  for (int i = 0; i <= n; i++) {
    for (int j = 0; j <= n; j++) {
      for (int k = 0; k <= n; k++) {
        for (int l = 0; l <= n; l++) {
          h_params_33(i,j,k,l).cut_xst_0 = cut_xst_0_33[i][j][k][l];
          h_params_33(i,j,k,l).cut_xst_c = cut_xst_c_33[i][j][k][l];
          h_params_33(i,j,k,l).cut_xst_lo = cut_xst_lo_33[i][j][k][l];
          h_params_33(i,j,k,l).cut_xst_hi = cut_xst_hi_33[i][j][k][l];
          h_params_33(i,j,k,l).cut_xst_lc = cut_xst_lc_33[i][j][k][l];
          h_params_33(i,j,k,l).cut_xst_hc = cut_xst_hc_33[i][j][k][l];
          h_params_33(i,j,k,l).a_xst4 = a_xst4_33[i][j][k][l];
          h_params_33(i,j,k,l).theta_xst4_0 = theta_xst4_0_33[i][j][k][l];
          h_params_33(i,j,k,l).dtheta_xst4_ast = dtheta_xst4_ast_33[i][j][k][l];
          h_params_33(i,j,k,l).b_xst4 = b_xst4_33[i][j][k][l];
          h_params_33(i,j,k,l).dtheta_xst4_c = dtheta_xst4_c_33[i][j][k][l];

          h_params_55(i,j,k,l).cut_xst_0 = cut_xst_0_55[i][j][k][l];
          h_params_55(i,j,k,l).cut_xst_c = cut_xst_c_55[i][j][k][l];
          h_params_55(i,j,k,l).cut_xst_lo = cut_xst_lo_55[i][j][k][l];
          h_params_55(i,j,k,l).cut_xst_hi = cut_xst_hi_55[i][j][k][l];
          h_params_55(i,j,k,l).cut_xst_lc = cut_xst_lc_55[i][j][k][l];
          h_params_55(i,j,k,l).cut_xst_hc = cut_xst_hc_55[i][j][k][l];
          h_params_55(i,j,k,l).a_xst4 = a_xst4_55[i][j][k][l];
          h_params_55(i,j,k,l).theta_xst4_0 = theta_xst4_0_55[i][j][k][l];
          h_params_55(i,j,k,l).dtheta_xst4_ast = dtheta_xst4_ast_55[i][j][k][l];
          h_params_55(i,j,k,l).b_xst4 = b_xst4_55[i][j][k][l];
          h_params_55(i,j,k,l).dtheta_xst4_c = dtheta_xst4_c_55[i][j][k][l];
        }
      }
    }
  }

  k_params_xstk.template modify<LMPHostType>();
  k_params_33.template modify<LMPHostType>();
  k_params_55.template modify<LMPHostType>();
  k_params_t7.template modify<LMPHostType>();
  k_params_t8.template modify<LMPHostType>();

  k_params_xstk.template sync<DeviceType>();
  k_params_33.template sync<DeviceType>();
  k_params_55.template sync<DeviceType>();
  k_params_t7.template sync<DeviceType>();
  k_params_t8.template sync<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxdna3XstkKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  // The eatom and vatom arrays are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
  auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  if (EFLAG) {
    if (eflag_atom) {
      const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * epair);
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx*fx);
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely*fy);
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz*fz);
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx*fy);
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx*fz);
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely*fz);

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
      } else {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        ev.v[1] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        ev.v[2] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        ev.v[3] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        ev.v[4] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        ev.v[5] += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
          a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
          a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
          a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
          a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
          a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
        if (NEWTON_PAIR || j < nlocal) {
        a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
        }
      } else {
        a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v0);
        a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v1);
        a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v2);
        a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v3);
        a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v4);
        a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * v5);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairOxdna3XstkKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairOxdna3XstkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdna3XstkKokkos<LMPHostType>;
#endif
}
