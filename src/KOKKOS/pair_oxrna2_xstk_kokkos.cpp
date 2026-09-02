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

#include "pair_oxrna2_xstk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "fix_oxdna_lrf_kokkos.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxrna2XstkKokkos<DeviceType>::PairOxrna2XstkKokkos(LAMMPS *lmp) : PairOxrna2Xstk(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  fix_oxdna_lrfKK = nullptr;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxrna2XstkKokkos<DeviceType>::~PairOxrna2XstkKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2XstkKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->create_kokkos(k_eatom, eatom, maxeatom, "pair:eatom");
    d_eatom = k_eatom.template view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom, vatom);
    memoryKK->create_kokkos(k_vatom, vatom, maxvatom, "pair:vatom");
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

  // get the neighbor list and neighbors used in operator()
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_neighbors = k_list->d_neighbors;
  anum = list->inum;
  d_alist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;

  int need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (eflag_atom) {
    if (need_dup) {
      dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    } else {
      ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    }
  }
  if (vflag_atom) {
    if (need_dup) {
      dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterDuplicated>(d_vatom);
    } else {
      ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
    }
  }

  copymode = 1;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->k_nx.template view<DeviceType>();
  d_ny_xtrct = fix_oxdna_lrfKK->k_ny.template view<DeviceType>();
  d_nz_xtrct = fix_oxdna_lrfKK->k_nz.template view<DeviceType>();

  // loop over pair interaction neighbors of my atoms
  EV_FLOAT ev;

  const int dispatch_neigh =
      (neighflag == HALF) ? 0 :
      (neighflag == HALFTHREAD) ? 1 :
      (neighflag == FULL) ? 2 : -1;

  if (dispatch_neigh < 0) {
    error->all(FLERR, "Unsupported neighbor flag in pair oxrna2/xstk/kk");
  }

  const int dispatch_key = (evflag ? 8 : 0) | (newton_pair ? 4 : 0) | dispatch_neigh;

  switch (dispatch_key) {
    case 0:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALF,0,0>>(0,anum),*this);
      break;
    case 1:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALFTHREAD,0,0>>(0,anum),*this);
      break;
    case 2:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<FULL,0,0>>(0,anum),*this);
      break;
    case 4:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALF,1,0>>(0,anum),*this);
      break;
    case 5:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALFTHREAD,1,0>>(0,anum),*this);
      break;
    case 6:
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<FULL,1,0>>(0,anum),*this);
      break;
    case 8:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALF,0,1>>(0,anum),*this,ev);
      break;
    case 9:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALFTHREAD,0,1>>(0,anum),*this,ev);
      break;
    case 10:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<FULL,0,1>>(0,anum),*this,ev);
      break;
    case 12:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALF,1,1>>(0,anum),*this,ev);
      break;
    case 13:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<HALFTHREAD,1,1>>(0,anum),*this,ev);
      break;
    case 14:
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairOxrna2XstkCompute<FULL,1,1>>(0,anum),*this,ev);
      break;
    default:
      error->all(FLERR, "Internal dispatch error in pair oxrna2/xstk/kk");
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
      Kokkos::Experimental::contribute(d_eatom,dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom,dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (need_dup) {
    if (eflag_atom) dup_eatom = decltype(dup_eatom)();
    if (vflag_atom) dup_vatom = decltype(dup_vatom)();
  } else {
    if (eflag_atom) ndup_eatom = decltype(ndup_eatom)();
    if (vflag_atom) ndup_vatom = decltype(ndup_vatom)();
  }

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxrna2XstkKokkos<DeviceType>::operator()(TagPairOxrna2XstkCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>,
  const int &ia, EV_FLOAT &ev) const
{
  Kokkos::View<KK_ACC_FLOAT *[3], typename AT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      a_f = f;
  Kokkos::View<KK_ACC_FLOAT *[3], typename AT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      a_torque = torque;

  const int a = d_alist(ia);
  const int atype = type(a);

  const KK_FLOAT a_nx0 = d_nx_xtrct(a, 0);
  const KK_FLOAT a_nx1 = d_nx_xtrct(a, 1);
  const KK_FLOAT a_nx2 = d_nx_xtrct(a, 2);
  const KK_FLOAT a_nz0 = d_nz_xtrct(a, 0);
  const KK_FLOAT a_nz1 = d_nz_xtrct(a, 1);
  const KK_FLOAT a_nz2 = d_nz_xtrct(a, 2);

  constexpr KK_FLOAT d_cbs = 0.4;
  const KK_FLOAT ra_cbs[3] = {d_cbs * a_nx0, d_cbs * a_nx1, d_cbs * a_nx2};

  const int bnum = d_numneigh(a);
  for (int ib = 0; ib < bnum; ib++) {
    int b = d_neighbors(a, ib);
    const KK_FLOAT factor_lj = special_lj[sbmask(b)];
    if (factor_lj == static_cast<KK_FLOAT>(0.0)) continue;
    b &= NEIGHMASK;

    const int btype = type(b);

    const KK_FLOAT b_nx0 = d_nx_xtrct(b, 0);
    const KK_FLOAT b_nx1 = d_nx_xtrct(b, 1);
    const KK_FLOAT b_nx2 = d_nx_xtrct(b, 2);
    const KK_FLOAT b_nz0 = d_nz_xtrct(b, 0);
    const KK_FLOAT b_nz1 = d_nz_xtrct(b, 1);
    const KK_FLOAT b_nz2 = d_nz_xtrct(b, 2);

    const KK_FLOAT rb_cbs[3] = {d_cbs * b_nx0, d_cbs * b_nx1, d_cbs * b_nx2};

    KK_FLOAT delr_bsbs[3];
    delr_bsbs[0] = x(a, 0) + ra_cbs[0] - x(b, 0) - rb_cbs[0];
    delr_bsbs[1] = x(a, 1) + ra_cbs[1] - x(b, 1) - rb_cbs[1];
    delr_bsbs[2] = x(a, 2) + ra_cbs[2] - x(b, 2) - rb_cbs[2];

    const KK_FLOAT rsq_bsbs = delr_bsbs[0] * delr_bsbs[0] + delr_bsbs[1] * delr_bsbs[1] +
        delr_bsbs[2] * delr_bsbs[2];
    if (rsq_bsbs <= static_cast<KK_FLOAT>(0.0)) continue;

    const KK_FLOAT r_bsbs = Kokkos::sqrt(rsq_bsbs);
    const KK_FLOAT rinv_bsbs = static_cast<KK_FLOAT>(1.0) / r_bsbs;

    KK_FLOAT delr_bsbs_norm[3];
    delr_bsbs_norm[0] = delr_bsbs[0] * rinv_bsbs;
    delr_bsbs_norm[1] = delr_bsbs[1] * rinv_bsbs;
    delr_bsbs_norm[2] = delr_bsbs[2] * rinv_bsbs;

    const KK_FLOAT f2 =
        F2_KK(r_bsbs, d_k_xst(atype, btype), d_cut_xst_0(atype, btype),
              d_cut_xst_lc(atype, btype), d_cut_xst_hc(atype, btype),
              d_cut_xst_lo(atype, btype), d_cut_xst_hi(atype, btype),
              d_b_xst_lo(atype, btype), d_b_xst_hi(atype, btype),
              d_cut_xst_c(atype, btype));
    if (f2 == static_cast<KK_FLOAT>(0.0)) continue;

    KK_FLOAT cost1 = -(a_nx0 * b_nx0 + a_nx1 * b_nx1 + a_nx2 * b_nx2);
    if (cost1 > 1.0) cost1 = 1.0;
    if (cost1 < -1.0) cost1 = -1.0;
    const KK_FLOAT theta1 = acos(cost1);
    const KK_FLOAT f4t1 = F4_KK(theta1, d_a_xst1(atype, btype), d_theta_xst1_0(atype, btype),
                                 d_dtheta_xst1_ast(atype, btype), d_b_xst1(atype, btype),
                                 d_dtheta_xst1_c(atype, btype));
    if (f4t1 == static_cast<KK_FLOAT>(0.0)) continue;

    KK_FLOAT cost2 = -(a_nx0 * delr_bsbs_norm[0] + a_nx1 * delr_bsbs_norm[1] +
                       a_nx2 * delr_bsbs_norm[2]);
    if (cost2 > 1.0) cost2 = 1.0;
    if (cost2 < -1.0) cost2 = -1.0;
    const KK_FLOAT theta2 = acos(cost2);
    const KK_FLOAT f4t2 = F4_KK(theta2, d_a_xst2(atype, btype), d_theta_xst2_0(atype, btype),
                                 d_dtheta_xst2_ast(atype, btype), d_b_xst2(atype, btype),
                                 d_dtheta_xst2_c(atype, btype));
    if (f4t2 == static_cast<KK_FLOAT>(0.0)) continue;

    KK_FLOAT cost3 = b_nx0 * delr_bsbs_norm[0] + b_nx1 * delr_bsbs_norm[1] +
        b_nx2 * delr_bsbs_norm[2];
    if (cost3 > 1.0) cost3 = 1.0;
    if (cost3 < -1.0) cost3 = -1.0;
    const KK_FLOAT theta3 = acos(cost3);
    const KK_FLOAT f4t3 = F4_KK(theta3, d_a_xst3(atype, btype), d_theta_xst3_0(atype, btype),
                                 d_dtheta_xst3_ast(atype, btype), d_b_xst3(atype, btype),
                                 d_dtheta_xst3_c(atype, btype));
    if (f4t3 == static_cast<KK_FLOAT>(0.0)) continue;

    KK_FLOAT cost7 = -(a_nz0 * delr_bsbs_norm[0] + a_nz1 * delr_bsbs_norm[1] +
                       a_nz2 * delr_bsbs_norm[2]);
    if (cost7 > 1.0) cost7 = 1.0;
    if (cost7 < -1.0) cost7 = -1.0;
    const KK_FLOAT theta7 = acos(cost7);
    const KK_FLOAT theta7p = static_cast<KK_FLOAT>(MathConst::MY_PI) - theta7;
    const KK_FLOAT f4t7 = F4_KK(theta7, d_a_xst7(atype, btype), d_theta_xst7_0(atype, btype),
                                 d_dtheta_xst7_ast(atype, btype), d_b_xst7(atype, btype),
                                 d_dtheta_xst7_c(atype, btype)) +
        F4_KK(theta7p, d_a_xst7(atype, btype), d_theta_xst7_0(atype, btype),
              d_dtheta_xst7_ast(atype, btype), d_b_xst7(atype, btype),
              d_dtheta_xst7_c(atype, btype));
    if (f4t7 == static_cast<KK_FLOAT>(0.0)) continue;

    KK_FLOAT cost8 = b_nz0 * delr_bsbs_norm[0] + b_nz1 * delr_bsbs_norm[1] +
        b_nz2 * delr_bsbs_norm[2];
    if (cost8 > 1.0) cost8 = 1.0;
    if (cost8 < -1.0) cost8 = -1.0;
    const KK_FLOAT theta8 = acos(cost8);
    const KK_FLOAT theta8p = static_cast<KK_FLOAT>(MathConst::MY_PI) - theta8;
    const KK_FLOAT f4t8 = F4_KK(theta8, d_a_xst8(atype, btype), d_theta_xst8_0(atype, btype),
                                 d_dtheta_xst8_ast(atype, btype), d_b_xst8(atype, btype),
                                 d_dtheta_xst8_c(atype, btype)) +
        F4_KK(theta8p, d_a_xst8(atype, btype), d_theta_xst8_0(atype, btype),
              d_dtheta_xst8_ast(atype, btype), d_b_xst8(atype, btype),
              d_dtheta_xst8_c(atype, btype));

    const KK_ACC_FLOAT evdwl = f2 * f4t1 * f4t2 * f4t3 * f4t7 * f4t8 * factor_lj;
    if (evdwl == static_cast<KK_FLOAT>(0.0)) continue;

    const KK_FLOAT df2 =
        DF2_KK(r_bsbs, d_k_xst(atype, btype), d_cut_xst_0(atype, btype),
               d_cut_xst_lc(atype, btype), d_cut_xst_hc(atype, btype),
               d_cut_xst_lo(atype, btype), d_cut_xst_hi(atype, btype),
               d_b_xst_lo(atype, btype), d_b_xst_hi(atype, btype));

    const KK_FLOAT df4t1 =
        DF4_KK(theta1, d_a_xst1(atype, btype), d_theta_xst1_0(atype, btype),
               d_dtheta_xst1_ast(atype, btype), d_b_xst1(atype, btype),
               d_dtheta_xst1_c(atype, btype)) /
        sin(theta1);

    const KK_FLOAT df4t2 =
        DF4_KK(theta2, d_a_xst2(atype, btype), d_theta_xst2_0(atype, btype),
               d_dtheta_xst2_ast(atype, btype), d_b_xst2(atype, btype),
               d_dtheta_xst2_c(atype, btype)) /
        sin(theta2);

    const KK_FLOAT df4t3 =
        DF4_KK(theta3, d_a_xst3(atype, btype), d_theta_xst3_0(atype, btype),
               d_dtheta_xst3_ast(atype, btype), d_b_xst3(atype, btype),
               d_dtheta_xst3_c(atype, btype)) /
        sin(theta3);

    const KK_FLOAT df4t7 =
        (DF4_KK(theta7, d_a_xst7(atype, btype), d_theta_xst7_0(atype, btype),
                d_dtheta_xst7_ast(atype, btype), d_b_xst7(atype, btype),
                d_dtheta_xst7_c(atype, btype)) -
         DF4_KK(theta7p, d_a_xst7(atype, btype), d_theta_xst7_0(atype, btype),
                d_dtheta_xst7_ast(atype, btype), d_b_xst7(atype, btype),
                d_dtheta_xst7_c(atype, btype))) /
        sin(theta7);

    const KK_FLOAT df4t8 =
        (DF4_KK(theta8, d_a_xst8(atype, btype), d_theta_xst8_0(atype, btype),
                d_dtheta_xst8_ast(atype, btype), d_b_xst8(atype, btype),
                d_dtheta_xst8_c(atype, btype)) -
         DF4_KK(theta8p, d_a_xst8(atype, btype), d_theta_xst8_0(atype, btype),
                d_dtheta_xst8_ast(atype, btype), d_b_xst8(atype, btype),
                d_dtheta_xst8_c(atype, btype))) /
        sin(theta8);

    KK_ACC_FLOAT delf[3] = {0.0, 0.0, 0.0};
    KK_ACC_FLOAT delta[3] = {0.0, 0.0, 0.0};
    KK_ACC_FLOAT deltb[3] = {0.0, 0.0, 0.0};

    KK_ACC_FLOAT finc = -df2 * f4t1 * f4t2 * f4t3 * f4t7 * f4t8 * rinv_bsbs * factor_lj;
    delf[0] += delr_bsbs[0] * finc;
    delf[1] += delr_bsbs[1] * finc;
    delf[2] += delr_bsbs[2] * finc;

    if (theta2 != static_cast<KK_FLOAT>(0.0)) {
      finc = -f2 * f4t1 * df4t2 * f4t3 * f4t7 * f4t8 * rinv_bsbs * factor_lj;
      delf[0] += (delr_bsbs_norm[0] * cost2 + a_nx0) * finc;
      delf[1] += (delr_bsbs_norm[1] * cost2 + a_nx1) * finc;
      delf[2] += (delr_bsbs_norm[2] * cost2 + a_nx2) * finc;
    }

    if (theta3 != static_cast<KK_FLOAT>(0.0)) {
      finc = -f2 * f4t1 * f4t2 * df4t3 * f4t7 * f4t8 * rinv_bsbs * factor_lj;
      delf[0] += (delr_bsbs_norm[0] * cost3 - b_nx0) * finc;
      delf[1] += (delr_bsbs_norm[1] * cost3 - b_nx1) * finc;
      delf[2] += (delr_bsbs_norm[2] * cost3 - b_nx2) * finc;
    }

    if (theta7 != static_cast<KK_FLOAT>(0.0)) {
      finc = -f2 * f4t1 * f4t2 * f4t3 * df4t7 * f4t8 * rinv_bsbs * factor_lj;
      delf[0] += (delr_bsbs_norm[0] * cost7 + a_nz0) * finc;
      delf[1] += (delr_bsbs_norm[1] * cost7 + a_nz1) * finc;
      delf[2] += (delr_bsbs_norm[2] * cost7 + a_nz2) * finc;
    }

    if (theta8 != static_cast<KK_FLOAT>(0.0)) {
      finc = -f2 * f4t1 * f4t2 * f4t3 * f4t7 * df4t8 * rinv_bsbs * factor_lj;
      delf[0] += (delr_bsbs_norm[0] * cost8 - b_nz0) * finc;
      delf[1] += (delr_bsbs_norm[1] * cost8 - b_nz1) * finc;
      delf[2] += (delr_bsbs_norm[2] * cost8 - b_nz2) * finc;
    }

    a_f(a, 0) += delf[0];
    a_f(a, 1) += delf[1];
    a_f(a, 2) += delf[2];

    delta[0] = ra_cbs[1] * delf[2] - ra_cbs[2] * delf[1];
    delta[1] = ra_cbs[2] * delf[0] - ra_cbs[0] * delf[2];
    delta[2] = ra_cbs[0] * delf[1] - ra_cbs[1] * delf[0];

    a_torque(a, 0) += delta[0];
    a_torque(a, 1) += delta[1];
    a_torque(a, 2) += delta[2];

    if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
      a_f(b, 0) -= delf[0];
      a_f(b, 1) -= delf[1];
      a_f(b, 2) -= delf[2];

      deltb[0] = rb_cbs[1] * delf[2] - rb_cbs[2] * delf[1];
      deltb[1] = rb_cbs[2] * delf[0] - rb_cbs[0] * delf[2];
      deltb[2] = rb_cbs[0] * delf[1] - rb_cbs[1] * delf[0];

      a_torque(b, 0) -= deltb[0];
      a_torque(b, 1) -= deltb[1];
      a_torque(b, 2) -= deltb[2];
    }

    if (EVFLAG) {
      ev.evdwl +=
          (((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || (b < nlocal)))
               ? 1.0
               : 0.5) *
          evdwl;

      if (vflag_either || eflag_atom) {
        this->template ev_tally_xyz<NEIGHFLAG, NEWTON_PAIR>(
            ev, a, b, evdwl, delf[0], delf[1], delf[2], x(a, 0) - x(b, 0),
            x(a, 1) - x(b, 1), x(a, 2) - x(b, 2));
      }
    }

    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;
    deltb[0] = 0.0;
    deltb[1] = 0.0;
    deltb[2] = 0.0;

    KK_ACC_FLOAT tpair;

    if (theta1 != static_cast<KK_FLOAT>(0.0)) {
      tpair = -f2 * df4t1 * f4t2 * f4t3 * f4t7 * f4t8 * factor_lj;
      const KK_FLOAT t1dir0 = a_nx1 * b_nx2 - a_nx2 * b_nx1;
      const KK_FLOAT t1dir1 = a_nx2 * b_nx0 - a_nx0 * b_nx2;
      const KK_FLOAT t1dir2 = a_nx0 * b_nx1 - a_nx1 * b_nx0;
      delta[0] += t1dir0 * tpair;
      delta[1] += t1dir1 * tpair;
      delta[2] += t1dir2 * tpair;
      deltb[0] += t1dir0 * tpair;
      deltb[1] += t1dir1 * tpair;
      deltb[2] += t1dir2 * tpair;
    }

    if (theta2 != static_cast<KK_FLOAT>(0.0)) {
      tpair = -f2 * f4t1 * df4t2 * f4t3 * f4t7 * f4t8 * factor_lj;
      const KK_FLOAT t2dir0 = a_nx1 * delr_bsbs_norm[2] - a_nx2 * delr_bsbs_norm[1];
      const KK_FLOAT t2dir1 = a_nx2 * delr_bsbs_norm[0] - a_nx0 * delr_bsbs_norm[2];
      const KK_FLOAT t2dir2 = a_nx0 * delr_bsbs_norm[1] - a_nx1 * delr_bsbs_norm[0];
      delta[0] += t2dir0 * tpair;
      delta[1] += t2dir1 * tpair;
      delta[2] += t2dir2 * tpair;
    }

    if (theta3 != static_cast<KK_FLOAT>(0.0)) {
      tpair = -f2 * f4t1 * f4t2 * df4t3 * f4t7 * f4t8 * factor_lj;
      const KK_FLOAT t3dir0 = b_nx1 * delr_bsbs_norm[2] - b_nx2 * delr_bsbs_norm[1];
      const KK_FLOAT t3dir1 = b_nx2 * delr_bsbs_norm[0] - b_nx0 * delr_bsbs_norm[2];
      const KK_FLOAT t3dir2 = b_nx0 * delr_bsbs_norm[1] - b_nx1 * delr_bsbs_norm[0];
      deltb[0] += t3dir0 * tpair;
      deltb[1] += t3dir1 * tpair;
      deltb[2] += t3dir2 * tpair;
    }

    if (theta7 != static_cast<KK_FLOAT>(0.0)) {
      tpair = -f2 * f4t1 * f4t2 * f4t3 * df4t7 * f4t8 * factor_lj;
      const KK_FLOAT t7dir0 = a_nz1 * delr_bsbs_norm[2] - a_nz2 * delr_bsbs_norm[1];
      const KK_FLOAT t7dir1 = a_nz2 * delr_bsbs_norm[0] - a_nz0 * delr_bsbs_norm[2];
      const KK_FLOAT t7dir2 = a_nz0 * delr_bsbs_norm[1] - a_nz1 * delr_bsbs_norm[0];
      delta[0] += t7dir0 * tpair;
      delta[1] += t7dir1 * tpair;
      delta[2] += t7dir2 * tpair;
    }

    if (theta8 != static_cast<KK_FLOAT>(0.0)) {
      tpair = -f2 * f4t1 * f4t2 * f4t3 * f4t7 * df4t8 * factor_lj;
      const KK_FLOAT t8dir0 = b_nz1 * delr_bsbs_norm[2] - b_nz2 * delr_bsbs_norm[1];
      const KK_FLOAT t8dir1 = b_nz2 * delr_bsbs_norm[0] - b_nz0 * delr_bsbs_norm[2];
      const KK_FLOAT t8dir2 = b_nz0 * delr_bsbs_norm[1] - b_nz1 * delr_bsbs_norm[0];
      deltb[0] += t8dir0 * tpair;
      deltb[1] += t8dir1 * tpair;
      deltb[2] += t8dir2 * tpair;
    }

    a_torque(a, 0) += delta[0];
    a_torque(a, 1) += delta[1];
    a_torque(a, 2) += delta[2];

    if ((NEIGHFLAG == HALF || NEIGHFLAG == HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
      a_torque(b, 0) -= deltb[0];
      a_torque(b, 1) -= deltb[1];
      a_torque(b, 2) -= deltb[2];
    }
  }
}

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxrna2XstkKokkos<DeviceType>::operator()(TagPairOxrna2XstkCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>,
  const int &ia) const
{
  EV_FLOAT ev;
  this->template operator()<NEIGHFLAG,NEWTON_PAIR,EVFLAG>
    (TagPairOxrna2XstkCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ia,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2XstkKokkos<DeviceType>::allocate()
{
  PairOxrna2Xstk::allocate();

  int n = atom->ntypes;

  memoryKK->create_kokkos(k_k_xst, n + 1, n + 1, "PairOxrna2Xstk:xst");
  memoryKK->create_kokkos(k_cut_xst_0, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_0");
  memoryKK->create_kokkos(k_cut_xst_c, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_c");
  memoryKK->create_kokkos(k_cut_xst_lo, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_lo");
  memoryKK->create_kokkos(k_cut_xst_hi, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_hi");
  memoryKK->create_kokkos(k_cut_xst_lc, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_lc");
  memoryKK->create_kokkos(k_cut_xst_hc, n + 1, n + 1, "PairOxrna2Xstk:cut_xst_hc");
  memoryKK->create_kokkos(k_b_xst_lo, n + 1, n + 1, "PairOxrna2Xstk:b_xst_lo");
  memoryKK->create_kokkos(k_b_xst_hi, n + 1, n + 1, "PairOxrna2Xstk:b_xst_hi");
  memoryKK->create_kokkos(k_cutsq_xst_hc, n + 1, n + 1, "PairOxrna2Xstk:cutsq_xst_hc");

  memoryKK->create_kokkos(k_a_xst1, n + 1, n + 1, "PairOxrna2Xstk:a_xst1");
  memoryKK->create_kokkos(k_theta_xst1_0, n + 1, n + 1, "PairOxrna2Xstk:theta_xst1_0");
  memoryKK->create_kokkos(k_dtheta_xst1_ast, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst1_ast");
  memoryKK->create_kokkos(k_b_xst1, n + 1, n + 1, "PairOxrna2Xstk:b_xst1");
  memoryKK->create_kokkos(k_dtheta_xst1_c, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst1_c");

  memoryKK->create_kokkos(k_a_xst2, n + 1, n + 1, "PairOxrna2Xstk:a_xst2");
  memoryKK->create_kokkos(k_theta_xst2_0, n + 1, n + 1, "PairOxrna2Xstk:theta_xst2_0");
  memoryKK->create_kokkos(k_dtheta_xst2_ast, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst2_ast");
  memoryKK->create_kokkos(k_b_xst2, n + 1, n + 1, "PairOxrna2Xstk:b_xst2");
  memoryKK->create_kokkos(k_dtheta_xst2_c, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst2_c");

  memoryKK->create_kokkos(k_a_xst3, n + 1, n + 1, "PairOxrna2Xstk:a_xst3");
  memoryKK->create_kokkos(k_theta_xst3_0, n + 1, n + 1, "PairOxrna2Xstk:theta_xst3_0");
  memoryKK->create_kokkos(k_dtheta_xst3_ast, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst3_ast");
  memoryKK->create_kokkos(k_b_xst3, n + 1, n + 1, "PairOxrna2Xstk:b_xst3");
  memoryKK->create_kokkos(k_dtheta_xst3_c, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst3_c");

  memoryKK->create_kokkos(k_a_xst7, n + 1, n + 1, "PairOxrna2Xstk:a_xst7");
  memoryKK->create_kokkos(k_theta_xst7_0, n + 1, n + 1, "PairOxrna2Xstk:theta_xst7_0");
  memoryKK->create_kokkos(k_dtheta_xst7_ast, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst7_ast");
  memoryKK->create_kokkos(k_b_xst7, n + 1, n + 1, "PairOxrna2Xstk:b_xst7");
  memoryKK->create_kokkos(k_dtheta_xst7_c, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst7_c");

  memoryKK->create_kokkos(k_a_xst8, n + 1, n + 1, "PairOxrna2Xstk:a_xst8");
  memoryKK->create_kokkos(k_theta_xst8_0, n + 1, n + 1, "PairOxrna2Xstk:theta_xst8_0");
  memoryKK->create_kokkos(k_dtheta_xst8_ast, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst8_ast");
  memoryKK->create_kokkos(k_b_xst8, n + 1, n + 1, "PairOxrna2Xstk:b_xst8");
  memoryKK->create_kokkos(k_dtheta_xst8_c, n + 1, n + 1, "PairOxrna2Xstk:dtheta_xst8_c");

  d_k_xst = k_k_xst.template view<DeviceType>();
  d_cut_xst_0 = k_cut_xst_0.template view<DeviceType>();
  d_cut_xst_c = k_cut_xst_c.template view<DeviceType>();
  d_cut_xst_lo = k_cut_xst_lo.template view<DeviceType>();
  d_cut_xst_hi = k_cut_xst_hi.template view<DeviceType>();
  d_cut_xst_lc = k_cut_xst_lc.template view<DeviceType>();
  d_cut_xst_hc = k_cut_xst_hc.template view<DeviceType>();
  d_b_xst_lo = k_b_xst_lo.template view<DeviceType>();
  d_b_xst_hi = k_b_xst_hi.template view<DeviceType>();
  d_cutsq_xst_hc = k_cutsq_xst_hc.template view<DeviceType>();

  d_a_xst1 = k_a_xst1.template view<DeviceType>();
  d_theta_xst1_0 = k_theta_xst1_0.template view<DeviceType>();
  d_dtheta_xst1_ast = k_dtheta_xst1_ast.template view<DeviceType>();
  d_b_xst1 = k_b_xst1.template view<DeviceType>();
  d_dtheta_xst1_c = k_dtheta_xst1_c.template view<DeviceType>();

  d_a_xst2 = k_a_xst2.template view<DeviceType>();
  d_theta_xst2_0 = k_theta_xst2_0.template view<DeviceType>();
  d_dtheta_xst2_ast = k_dtheta_xst2_ast.template view<DeviceType>();
  d_b_xst2 = k_b_xst2.template view<DeviceType>();
  d_dtheta_xst2_c = k_dtheta_xst2_c.template view<DeviceType>();

  d_a_xst3 = k_a_xst3.template view<DeviceType>();
  d_theta_xst3_0 = k_theta_xst3_0.template view<DeviceType>();
  d_dtheta_xst3_ast = k_dtheta_xst3_ast.template view<DeviceType>();
  d_b_xst3 = k_b_xst3.template view<DeviceType>();
  d_dtheta_xst3_c = k_dtheta_xst3_c.template view<DeviceType>();

  d_a_xst7 = k_a_xst7.template view<DeviceType>();
  d_theta_xst7_0 = k_theta_xst7_0.template view<DeviceType>();
  d_dtheta_xst7_ast = k_dtheta_xst7_ast.template view<DeviceType>();
  d_b_xst7 = k_b_xst7.template view<DeviceType>();
  d_dtheta_xst7_c = k_dtheta_xst7_c.template view<DeviceType>();

  d_a_xst8 = k_a_xst8.template view<DeviceType>();
  d_theta_xst8_0 = k_theta_xst8_0.template view<DeviceType>();
  d_dtheta_xst8_ast = k_dtheta_xst8_ast.template view<DeviceType>();
  d_b_xst8 = k_b_xst8.template view<DeviceType>();
  d_dtheta_xst8_c = k_dtheta_xst8_c.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2XstkKokkos<DeviceType>::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2XstkKokkos<DeviceType>::init_style()
{
  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0)
    error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else
    fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);
  if (!fix_oxdna_lrfKK)
    error->all(FLERR, "Fix OXDNA/LRF/kk lookup failed");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxrna2XstkKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxrna2Xstk::init_one(i, j);

  k_k_xst.view_host()(i, j) = k_xst[i][j];
  k_k_xst.view_host()(j, i) = k_xst[j][i];
  k_cut_xst_0.view_host()(i, j) = cut_xst_0[i][j];
  k_cut_xst_0.view_host()(j, i) = cut_xst_0[j][i];
  k_cut_xst_c.view_host()(i, j) = cut_xst_c[i][j];
  k_cut_xst_c.view_host()(j, i) = cut_xst_c[j][i];
  k_cut_xst_lo.view_host()(i, j) = cut_xst_lo[i][j];
  k_cut_xst_lo.view_host()(j, i) = cut_xst_lo[j][i];
  k_cut_xst_hi.view_host()(i, j) = cut_xst_hi[i][j];
  k_cut_xst_hi.view_host()(j, i) = cut_xst_hi[j][i];
  k_cut_xst_lc.view_host()(i, j) = cut_xst_lc[i][j];
  k_cut_xst_lc.view_host()(j, i) = cut_xst_lc[j][i];
  k_cut_xst_hc.view_host()(i, j) = cut_xst_hc[i][j];
  k_cut_xst_hc.view_host()(j, i) = cut_xst_hc[j][i];
  k_b_xst_lo.view_host()(i, j) = b_xst_lo[i][j];
  k_b_xst_lo.view_host()(j, i) = b_xst_lo[j][i];
  k_b_xst_hi.view_host()(i, j) = b_xst_hi[i][j];
  k_b_xst_hi.view_host()(j, i) = b_xst_hi[j][i];
  k_cutsq_xst_hc.view_host()(i, j) = cutsq_xst_hc[i][j];
  k_cutsq_xst_hc.view_host()(j, i) = cutsq_xst_hc[j][i];

  k_a_xst1.view_host()(i, j) = a_xst1[i][j];
  k_a_xst1.view_host()(j, i) = a_xst1[j][i];
  k_theta_xst1_0.view_host()(i, j) = theta_xst1_0[i][j];
  k_theta_xst1_0.view_host()(j, i) = theta_xst1_0[j][i];
  k_dtheta_xst1_ast.view_host()(i, j) = dtheta_xst1_ast[i][j];
  k_dtheta_xst1_ast.view_host()(j, i) = dtheta_xst1_ast[j][i];
  k_b_xst1.view_host()(i, j) = b_xst1[i][j];
  k_b_xst1.view_host()(j, i) = b_xst1[j][i];
  k_dtheta_xst1_c.view_host()(i, j) = dtheta_xst1_c[i][j];
  k_dtheta_xst1_c.view_host()(j, i) = dtheta_xst1_c[j][i];

  k_a_xst2.view_host()(i, j) = a_xst2[i][j];
  k_a_xst2.view_host()(j, i) = a_xst2[j][i];
  k_theta_xst2_0.view_host()(i, j) = theta_xst2_0[i][j];
  k_theta_xst2_0.view_host()(j, i) = theta_xst2_0[j][i];
  k_dtheta_xst2_ast.view_host()(i, j) = dtheta_xst2_ast[i][j];
  k_dtheta_xst2_ast.view_host()(j, i) = dtheta_xst2_ast[j][i];
  k_b_xst2.view_host()(i, j) = b_xst2[i][j];
  k_b_xst2.view_host()(j, i) = b_xst2[j][i];
  k_dtheta_xst2_c.view_host()(i, j) = dtheta_xst2_c[i][j];
  k_dtheta_xst2_c.view_host()(j, i) = dtheta_xst2_c[j][i];

  k_a_xst3.view_host()(i, j) = a_xst3[i][j];
  k_a_xst3.view_host()(j, i) = a_xst3[j][i];
  k_theta_xst3_0.view_host()(i, j) = theta_xst3_0[i][j];
  k_theta_xst3_0.view_host()(j, i) = theta_xst3_0[j][i];
  k_dtheta_xst3_ast.view_host()(i, j) = dtheta_xst3_ast[i][j];
  k_dtheta_xst3_ast.view_host()(j, i) = dtheta_xst3_ast[j][i];
  k_b_xst3.view_host()(i, j) = b_xst3[i][j];
  k_b_xst3.view_host()(j, i) = b_xst3[j][i];
  k_dtheta_xst3_c.view_host()(i, j) = dtheta_xst3_c[i][j];
  k_dtheta_xst3_c.view_host()(j, i) = dtheta_xst3_c[j][i];

  k_a_xst7.view_host()(i, j) = a_xst7[i][j];
  k_a_xst7.view_host()(j, i) = a_xst7[j][i];
  k_theta_xst7_0.view_host()(i, j) = theta_xst7_0[i][j];
  k_theta_xst7_0.view_host()(j, i) = theta_xst7_0[j][i];
  k_dtheta_xst7_ast.view_host()(i, j) = dtheta_xst7_ast[i][j];
  k_dtheta_xst7_ast.view_host()(j, i) = dtheta_xst7_ast[j][i];
  k_b_xst7.view_host()(i, j) = b_xst7[i][j];
  k_b_xst7.view_host()(j, i) = b_xst7[j][i];
  k_dtheta_xst7_c.view_host()(i, j) = dtheta_xst7_c[i][j];
  k_dtheta_xst7_c.view_host()(j, i) = dtheta_xst7_c[j][i];

  k_a_xst8.view_host()(i, j) = a_xst8[i][j];
  k_a_xst8.view_host()(j, i) = a_xst8[j][i];
  k_theta_xst8_0.view_host()(i, j) = theta_xst8_0[i][j];
  k_theta_xst8_0.view_host()(j, i) = theta_xst8_0[j][i];
  k_dtheta_xst8_ast.view_host()(i, j) = dtheta_xst8_ast[i][j];
  k_dtheta_xst8_ast.view_host()(j, i) = dtheta_xst8_ast[j][i];
  k_b_xst8.view_host()(i, j) = b_xst8[i][j];
  k_b_xst8.view_host()(j, i) = b_xst8[j][i];
  k_dtheta_xst8_c.view_host()(i, j) = dtheta_xst8_c[i][j];
  k_dtheta_xst8_c.view_host()(j, i) = dtheta_xst8_c[j][i];

  k_k_xst.template modify<LMPHostType>();
  k_cut_xst_0.template modify<LMPHostType>();
  k_cut_xst_c.template modify<LMPHostType>();
  k_cut_xst_lo.template modify<LMPHostType>();
  k_cut_xst_hi.template modify<LMPHostType>();
  k_cut_xst_lc.template modify<LMPHostType>();
  k_cut_xst_hc.template modify<LMPHostType>();
  k_b_xst_lo.template modify<LMPHostType>();
  k_b_xst_hi.template modify<LMPHostType>();
  k_cutsq_xst_hc.template modify<LMPHostType>();

  k_a_xst1.template modify<LMPHostType>();
  k_theta_xst1_0.template modify<LMPHostType>();
  k_dtheta_xst1_ast.template modify<LMPHostType>();
  k_b_xst1.template modify<LMPHostType>();
  k_dtheta_xst1_c.template modify<LMPHostType>();

  k_a_xst2.template modify<LMPHostType>();
  k_theta_xst2_0.template modify<LMPHostType>();
  k_dtheta_xst2_ast.template modify<LMPHostType>();
  k_b_xst2.template modify<LMPHostType>();
  k_dtheta_xst2_c.template modify<LMPHostType>();

  k_a_xst3.template modify<LMPHostType>();
  k_theta_xst3_0.template modify<LMPHostType>();
  k_dtheta_xst3_ast.template modify<LMPHostType>();
  k_b_xst3.template modify<LMPHostType>();
  k_dtheta_xst3_c.template modify<LMPHostType>();

  k_a_xst7.template modify<LMPHostType>();
  k_theta_xst7_0.template modify<LMPHostType>();
  k_dtheta_xst7_ast.template modify<LMPHostType>();
  k_b_xst7.template modify<LMPHostType>();
  k_dtheta_xst7_c.template modify<LMPHostType>();

  k_a_xst8.template modify<LMPHostType>();
  k_theta_xst8_0.template modify<LMPHostType>();
  k_dtheta_xst8_ast.template modify<LMPHostType>();
  k_b_xst8.template modify<LMPHostType>();
  k_dtheta_xst8_c.template modify<LMPHostType>();

  k_k_xst.template sync<DeviceType>();
  k_cut_xst_0.template sync<DeviceType>();
  k_cut_xst_c.template sync<DeviceType>();
  k_cut_xst_lo.template sync<DeviceType>();
  k_cut_xst_hi.template sync<DeviceType>();
  k_cut_xst_lc.template sync<DeviceType>();
  k_cut_xst_hc.template sync<DeviceType>();
  k_b_xst_lo.template sync<DeviceType>();
  k_b_xst_hi.template sync<DeviceType>();
  k_cutsq_xst_hc.template sync<DeviceType>();

  k_a_xst1.template sync<DeviceType>();
  k_theta_xst1_0.template sync<DeviceType>();
  k_dtheta_xst1_ast.template sync<DeviceType>();
  k_b_xst1.template sync<DeviceType>();
  k_dtheta_xst1_c.template sync<DeviceType>();

  k_a_xst2.template sync<DeviceType>();
  k_theta_xst2_0.template sync<DeviceType>();
  k_dtheta_xst2_ast.template sync<DeviceType>();
  k_b_xst2.template sync<DeviceType>();
  k_dtheta_xst2_c.template sync<DeviceType>();

  k_a_xst3.template sync<DeviceType>();
  k_theta_xst3_0.template sync<DeviceType>();
  k_dtheta_xst3_ast.template sync<DeviceType>();
  k_b_xst3.template sync<DeviceType>();
  k_dtheta_xst3_c.template sync<DeviceType>();

  k_a_xst7.template sync<DeviceType>();
  k_theta_xst7_0.template sync<DeviceType>();
  k_dtheta_xst7_ast.template sync<DeviceType>();
  k_b_xst7.template sync<DeviceType>();
  k_dtheta_xst7_c.template sync<DeviceType>();

  k_a_xst8.template sync<DeviceType>();
  k_theta_xst8_0.template sync<DeviceType>();
  k_dtheta_xst8_ast.template sync<DeviceType>();
  k_b_xst8.template sync<DeviceType>();
  k_dtheta_xst8_c.template sync<DeviceType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxrna2XstkKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
  const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
  const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  const int EFLAG = eflag;
  const int VFLAG = vflag_either;

  if (EFLAG) {
    if (eflag_atom) {
      // The eatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial.
      auto v_eatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,
        decltype(dup_eatom),decltype(ndup_eatom)>::get(dup_eatom,ndup_eatom);
      auto a_eatom = v_eatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

      const KK_ACC_FLOAT epairhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * epair);
      if (NEIGHFLAG != FULL) {
        if (NEWTON_PAIR || i < nlocal) a_eatom[i] += epairhalf;
        if (NEWTON_PAIR || j < nlocal) a_eatom[j] += epairhalf;
      } else {
        a_eatom[i] += epairhalf;
      }
    }
  }

  if (VFLAG) {
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx * fx);
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely * fy);
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz * fz);
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx * fy);
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx * fz);
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely * fz);

    if (vflag_global) {
      if (NEIGHFLAG != FULL) {
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
      // The vatom array is duplicated for OpenMP, atomic for GPU, and neither for Serial.
      auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,
        decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
      auto a_vatom = v_vatom.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

      if (NEIGHFLAG != FULL) {
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
KOKKOS_INLINE_FUNCTION int PairOxrna2XstkKokkos<DeviceType>::sbmask(const int &j) const
{
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairOxrna2XstkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxrna2XstkKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
