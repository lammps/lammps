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

#include "pair_oxdna_hbond_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_npair_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;

// NOTE: I've introduced some extra early returns in calls related to ComputeGPUPair.
// With the use of fma and trig identity "sin^2(theta) = 1 - cos^2(theta)", some of the
// math ops yeild unstable/seg-fault results without these - especially when running
// FP32.

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaHbondKokkos<DeviceType>::PairOxdnaHbondKokkos(LAMMPS *lmp) : PairOxdnaHbond(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK | TAG_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA;
  screened_pair_count = 0;
  unique_basepair_enabled = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaHbondKokkos<DeviceType>::~PairOxdnaHbondKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  tag = atomKK->k_tag.template view<DeviceType>();

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

  // If we're on a GPU, look up fix_oxdna_npairKK screened pair count and pair_a/b views.
  if (execution_space != HostKK) {
    screened_pair_count = fix_oxdna_npairKK->screened_pair_count;
    d_pairs_screened = fix_oxdna_npairKK->k_pairs_screened.template view<DeviceType>();
  }

  if (unique_basepair_enabled && idc) {
    const int nall = atom->nlocal + atom->nghost;
    if (k_idc.extent(0) < static_cast<size_t>(nall)) {
      k_idc = DAT::tdual_int_1d("pair:idc", nall);
    }
    auto h_idc = k_idc.view_host();
    for (int i = 0; i < nall; i++) h_idc(i) = idc[i];
    k_idc.template modify<LMPHostType>();
    k_idc.template sync<DeviceType>();
    d_idc = k_idc.template view<DeviceType>();
  }

  // loop over neighbors of my atoms for compute functors

  EV_FLOAT ev;

  // Host/GPU launch paths are split to avoid an execution-space branch per dispatch call.
  auto run_compute_host = [&](auto host_tag, auto evflag_tag) {
    constexpr int EVFLAG = decltype(evflag_tag)::value;
    if constexpr (EVFLAG) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, decltype(host_tag)>(0,anum),*this,ev);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, decltype(host_tag)>(0,anum),*this);
    }
  };
  auto run_compute_gpu = [&](auto gpu_tag, auto evflag_tag) {
    constexpr int EVFLAG = decltype(evflag_tag)::value;
    if constexpr (EVFLAG) {
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_pair_count),*this,ev);
    } else {
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_pair_count),*this);
    }
  };

  const bool use_host_launch = (execution_space == HostKK);

  auto run_compute_by_flags = [&](auto neighflag_tag, auto newtonpair_tag, auto evflag_tag) {
    constexpr int NEIGHFLAG = decltype(neighflag_tag)::value;
    constexpr int NEWTON_PAIR = decltype(newtonpair_tag)::value;
    constexpr int EVFLAG = decltype(evflag_tag)::value;

    if (oxdnaflag == OXDNA) {
      if (use_host_launch) {
        run_compute_host(TagPairOxdnaHbondCompute<OXDNA,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      } else {
        run_compute_gpu(TagPairOxdnaHbondComputeGPUPair<OXDNA,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      }
    } else if (oxdnaflag == OXDNA3) {
      if (use_host_launch) {
        run_compute_host(TagPairOxdnaHbondCompute<OXDNA3,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      } else {
        run_compute_gpu(TagPairOxdnaHbondComputeGPUPair<OXDNA3,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      }
    } else {
      error->all(FLERR, "Unknown OXDNA model flag in pair oxdna/hbond/kk");
    }
  };

  const int dispatch_neigh =
      (neighflag == HALF) ? 0 :
      (neighflag == HALFTHREAD) ? 1 :
      (neighflag == FULL) ? 2 : -1;

  if (dispatch_neigh < 0) {
    error->all(FLERR, "Unsupported neighbor flag in pair oxdna/hbond/kk");
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
    default: error->all(FLERR, "Internal dispatch error in pair oxdna/hbond/kk");
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
   Standard non-GPU Compute Functor(s)
-------------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia, EV_FLOAT &ev) const
{
  // f and torque array are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  const int a = d_alist(ia);
  const int atype = type(a);
  // vectors COM-hbond site in lab frame
  KK_FLOAT ra_chb[3], rb_chb[3];

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];    // force, torque increment
  KK_ACC_FLOAT evdwl, finc, tpair;             // energy, force, torque
  KK_FLOAT delr_hb[3],delr_hb_norm[3],rsq_hb,r_hb,rinv_hb;
  KK_FLOAT theta1,t1dir[3],cost1;
  KK_FLOAT theta2,t2dir[3],cost2;
  KK_FLOAT theta3,t3dir[3],cost3;
  KK_FLOAT theta4,t4dir[3],cost4;
  KK_FLOAT theta7,t7dir[3],cost7;
  KK_FLOAT theta8,t8dir[3],cost8;

  KK_FLOAT f1,f4t1,f4t4,f4t2,f4t3,f4t7,f4t8;
  KK_FLOAT df1,df4t1,df4t4,df4t2,df4t3,df4t7,df4t8;

  // vector COM-hbond site a
  if (OXDNAFLAG == OXDNA) {
    // Used for oxdna1 and oxdna2, which have the same hbond site offset.
    constexpr KK_FLOAT dx_cbs_oxdna1=+0.4;
    ra_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
  } else if (OXDNAFLAG == OXDNA3) {
    constexpr KK_FLOAT dx_cbs_pur_oxdna3 = +0.43;
    constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = +0.37;
    int nucl_acid = (atype%4);
    if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
      ra_chb[0] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,0);
      ra_chb[1] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,1);
      ra_chb[2] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,2);
    } else {  // purine (A or G)
      ra_chb[0] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,0);
      ra_chb[1] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,1);
      ra_chb[2] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,2);
    }
  }

  const int bnum = d_numneigh(a);

  for (int ib = 0; ib < bnum; ib++) {

    int b = d_neighbors(a,ib);
    const KK_FLOAT factor_lj = special_lj[sbmask(b)];
    b &= NEIGHMASK;
    const int btype = type(b);

    if (unique_basepair_enabled) {
      const int idca = d_idc(a);
      const int idcb = d_idc(b);
      if (idca != tag(b) && idcb != tag(a) && idca > 0 && idcb > 0) continue;
    }

    // vector COM-hbond site b
    if (OXDNAFLAG == OXDNA) {
      // Used for oxdna1 and oxdna2, which have the same hbond site offset.
      constexpr KK_FLOAT dx_cbs_oxdna1=+0.4;
      rb_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(b,0);
      rb_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(b,1);
      rb_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(b,2);
    } else if (OXDNAFLAG == OXDNA3) {
      constexpr KK_FLOAT dx_cbs_pur_oxdna3 = +0.43;
      constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = +0.37;
      int nucl_acid = (btype%4);
      if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
        rb_chb[0] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,0);
        rb_chb[1] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,1);
        rb_chb[2] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,2);
      } else {  // purine (A or G)
        rb_chb[0] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,0);
        rb_chb[1] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,1);
        rb_chb[2] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,2);
      }
    }

    // vector h-bonding site b-a
    delr_hb[0] = x(a,0) + ra_chb[0] - x(b,0) - rb_chb[0];
    delr_hb[1] = x(a,1) + ra_chb[1] - x(b,1) - rb_chb[1];
    delr_hb[2] = x(a,2) + ra_chb[2] - x(b,2) - rb_chb[2];

    rsq_hb = delr_hb[0]*delr_hb[0] + delr_hb[1]*delr_hb[1] + delr_hb[2]*delr_hb[2];
    r_hb = sqrtf(rsq_hb);
    rinv_hb = 1.0 / r_hb;

    delr_hb_norm[0] = delr_hb[0] * rinv_hb;
    delr_hb_norm[1] = delr_hb[1] * rinv_hb;
    delr_hb_norm[2] = delr_hb[2] * rinv_hb;

    // beginning of modulation factors

    // f1 = f1 modulation factor
    f1 = F1_KK(r_hb, d_epsilon_hb(atype,btype), d_a_hb(atype,btype), d_cut_hb_0(atype,btype),
            d_cut_hb_lc(atype,btype), d_cut_hb_hc(atype,btype), d_cut_hb_lo(atype,btype),
            d_cut_hb_hi(atype,btype), d_b_hb_lo(atype,btype),
            d_b_hb_hi(atype,btype), d_shift_hb(atype,btype));

    // start early rejection criterium
    if (f1) {
      // theta1 calculation
      cost1 = - (d_nx_xtrct(a,0)*d_nx_xtrct(b,0) + d_nx_xtrct(a,1)*d_nx_xtrct(b,1) + d_nx_xtrct(a,2)*d_nx_xtrct(b,2));
      if (cost1 > 1.0) cost1 = 1.0;
      if (cost1 < -1.0) cost1 = -1.0;
      theta1 = acos(cost1);
      // f4t1 = f4 modulation factor
      f4t1 = F4_KK(theta1, d_a_hb1(atype,btype), d_theta_hb1_0(atype, btype), d_dtheta_hb1_ast(atype, btype),
              d_b_hb1(atype, btype), d_dtheta_hb1_c(atype, btype));
    // end of f1

    // f4t1 early rejection criterium
    if (f4t1) {
      // theta2 calculation
      cost2 = - (d_nx_xtrct(a,0)*delr_hb_norm[0] + d_nx_xtrct(a,1)*delr_hb_norm[1] + d_nx_xtrct(a,2)*delr_hb_norm[2]);
      if (cost2 > 1.0) cost2 = 1.0;
      if (cost2 < -1.0) cost2 = -1.0;
      theta2 = acos(cost2);
      // f4t2 = f4 modulation factor
      f4t2 = F4_KK(theta2, d_a_hb2(atype,btype), d_theta_hb2_0(atype, btype), d_dtheta_hb2_ast(atype, btype),
              d_b_hb2(atype, btype), d_dtheta_hb2_c(atype, btype));
    // end of f4t1

    // f4t2 early rejection criterium
    if (f4t2) {
      // theta3 calculation
      cost3 = d_nx_xtrct(b,0)*delr_hb_norm[0] + d_nx_xtrct(b,1)*delr_hb_norm[1] + d_nx_xtrct(b,2)*delr_hb_norm[2];
      if (cost3 > 1.0) cost3 = 1.0;
      if (cost3 < -1.0) cost3 = -1.0;
      theta3 = acos(cost3);
      // f4t3 = f4 modulation factor
      f4t3 = F4_KK(theta3, d_a_hb3(atype,btype), d_theta_hb3_0(atype, btype), d_dtheta_hb3_ast(atype, btype),
              d_b_hb3(atype, btype), d_dtheta_hb3_c(atype, btype));
    // end of f4t2

    // f4t3 early rejection criterium
    if (f4t3) {
      // theta4 calculation
      cost4 = d_nz_xtrct(a,0)*d_nz_xtrct(b,0) + d_nz_xtrct(a,1)*d_nz_xtrct(b,1) + d_nz_xtrct(a,2)*d_nz_xtrct(b,2);
      if (cost4 > 1.0) cost4 = 1.0;
      if (cost4 < -1.0) cost4 = -1.0;
      theta4 = acos(cost4);
      // f4t4 = f4 modulation factor
      f4t4 = F4_KK(theta4, d_a_hb4(atype,btype), d_theta_hb4_0(atype, btype), d_dtheta_hb4_ast(atype, btype),
              d_b_hb4(atype, btype), d_dtheta_hb4_c(atype, btype));
    // end of f4t3

    // f4t4 early rejection criterium
    if (f4t4) {
      cost7 = - (d_nz_xtrct(a,0)*delr_hb_norm[0] + d_nz_xtrct(a,1)*delr_hb_norm[1] + d_nz_xtrct(a,2)*delr_hb_norm[2]);
      if (cost7 > 1.0) cost7 = 1.0;
      if (cost7 < -1.0) cost7 = -1.0;
      theta7 = acos(cost7);
      // f4t7 = f4 modulation factor
      f4t7 = F4_KK(theta7, d_a_hb7(atype,btype), d_theta_hb7_0(atype, btype), d_dtheta_hb7_ast(atype, btype),
              d_b_hb7(atype, btype), d_dtheta_hb7_c(atype, btype));
    // end of f4t4

    // f4t7 early rejection criterium
    if (f4t7) {
      cost8 = d_nz_xtrct(b,0)*delr_hb_norm[0] + d_nz_xtrct(b,1)*delr_hb_norm[1] + d_nz_xtrct(b,2)*delr_hb_norm[2];
      if (cost8 > 1.0) cost8 = 1.0;
      if (cost8 < -1.0) cost8 = -1.0;
      theta8 = acos(cost8);
      // f4t8 = f4 modulation factor
      f4t8 = F4_KK(theta8, d_a_hb8(atype,btype), d_theta_hb8_0(atype, btype), d_dtheta_hb8_ast(atype, btype),
              d_b_hb8(atype, btype), d_dtheta_hb8_c(atype, btype));

      evdwl = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
    // end of f4t7

    // evdwl early rejection criterium
    if (evdwl) {
      // df1 = DF1 modulation factor
      df1 = DF1_KK(r_hb, d_epsilon_hb(atype,btype), d_a_hb(atype,btype), d_cut_hb_0(atype,btype),
            d_cut_hb_lc(atype,btype), d_cut_hb_hc(atype,btype), d_cut_hb_lo(atype,btype),
            d_cut_hb_hi(atype,btype), d_b_hb_lo(atype,btype),
            d_b_hb_hi(atype,btype));
      // df4t1 = DF4 modulation factor
      df4t1 = DF4_KK(theta1, d_a_hb1(atype,btype), d_theta_hb1_0(atype, btype), d_dtheta_hb1_ast(atype, btype),
              d_b_hb1(atype, btype), d_dtheta_hb1_c(atype, btype))/sin(theta1);
      // df4t2 = DF4 modulation factor
      df4t2 = DF4_KK(theta2, d_a_hb2(atype,btype), d_theta_hb2_0(atype, btype), d_dtheta_hb2_ast(atype, btype),
              d_b_hb2(atype, btype), d_dtheta_hb2_c(atype, btype))/sin(theta2);
      // df4t3 = DF4 modulation factor
      df4t3 = DF4_KK(theta3, d_a_hb3(atype,btype), d_theta_hb3_0(atype, btype), d_dtheta_hb3_ast(atype, btype),
              d_b_hb3(atype, btype), d_dtheta_hb3_c(atype, btype))/sin(theta3);
      // df4t4 = DF4 modulation factor
      df4t4 = DF4_KK(theta4, d_a_hb4(atype,btype), d_theta_hb4_0(atype, btype), d_dtheta_hb4_ast(atype, btype),
              d_b_hb4(atype, btype), d_dtheta_hb4_c(atype, btype))/sin(theta4);
      // df4t7 = DF4 modulation factor
      df4t7 = DF4_KK(theta7, d_a_hb7(atype,btype), d_theta_hb7_0(atype, btype), d_dtheta_hb7_ast(atype, btype),
              d_b_hb7(atype, btype), d_dtheta_hb7_c(atype, btype))/sin(theta7);
      // df4t8 = DF4 modulation factor
      df4t8 = DF4_KK(theta8, d_a_hb8(atype,btype), d_theta_hb8_0(atype, btype), d_dtheta_hb8_ast(atype, btype),
              d_b_hb8(atype, btype), d_dtheta_hb8_c(atype, btype))/sin(theta8);

      // force, torque, and viral contributions for forces between h-bonding sites

      delf[0] = 0.0;
      delf[1] = 0.0;
      delf[2] = 0.0;

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;

      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // radial force
      finc  = -df1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

      delf[0] += delr_hb[0] * finc;
      delf[1] += delr_hb[1] * finc;
      delf[2] += delr_hb[2] * finc;

      // theta2 force
      if (theta2 != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;

        delf[0] += (delr_hb_norm[0]*cost2 + d_nx_xtrct(a,0)) * finc;
        delf[1] += (delr_hb_norm[1]*cost2 + d_nx_xtrct(a,1)) * finc;
        delf[2] += (delr_hb_norm[2]*cost2 + d_nx_xtrct(a,2)) * finc;
      }

      // theta3 force
      if (theta3 != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;

        delf[0] += (delr_hb_norm[0]*cost3 - d_nx_xtrct(b,0)) * finc;
        delf[1] += (delr_hb_norm[1]*cost3 - d_nx_xtrct(b,1)) * finc;
        delf[2] += (delr_hb_norm[2]*cost3 - d_nx_xtrct(b,2)) * finc;
      }

      // theta7 force
      if (theta7 != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * rinv_hb * factor_lj;

        delf[0] += (delr_hb_norm[0]*cost7 + d_nz_xtrct(a,0)) * finc;
        delf[1] += (delr_hb_norm[1]*cost7 + d_nz_xtrct(a,1)) * finc;
        delf[2] += (delr_hb_norm[2]*cost7 + d_nz_xtrct(a,2)) * finc;

      }

      // theta8 force
      if (theta8 != static_cast<KK_FLOAT>(0.0)) {

        finc  = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * rinv_hb * factor_lj;

        delf[0] += (delr_hb_norm[0]*cost8 - d_nz_xtrct(b,0)) * finc;
        delf[1] += (delr_hb_norm[1]*cost8 - d_nz_xtrct(b,1)) * finc;
        delf[2] += (delr_hb_norm[2]*cost8 - d_nz_xtrct(b,2)) * finc;

      }

      // increment forces and torques

      a_f(a,0) += delf[0];
      a_f(a,1) += delf[1];
      a_f(a,2) += delf[2];
      delta[0] = ra_chb[1]*delf[2] - ra_chb[2]*delf[1];
      delta[1] = ra_chb[2]*delf[0] - ra_chb[0]*delf[2];
      delta[2] = ra_chb[0]*delf[1] - ra_chb[1]*delf[0];
      a_torque(a,0) += delta[0];
      a_torque(a,1) += delta[1];
      a_torque(a,2) += delta[2];

      if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
        a_f(b,0) -= delf[0];
        a_f(b,1) -= delf[1];
        a_f(b,2) -= delf[2];
        deltb[0] = rb_chb[1]*delf[2] - rb_chb[2]*delf[1];
        deltb[1] = rb_chb[2]*delf[0] - rb_chb[0]*delf[2];
        deltb[2] = rb_chb[0]*delf[1] - rb_chb[1]*delf[0];
        a_torque(b,0) -= deltb[0];
        a_torque(b,1) -= deltb[1];
        a_torque(b,2) -= deltb[2];
      }

      // increment energy and virial
      // NOTE: The virial is calculated on the 'molecular' basis.
      // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

      if (EVFLAG) {
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl,\
          delf[0],delf[1],delf[2],x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
        }
      }

      // pure torques not expressible as r x f

      delta[0] = 0.0;
      delta[1] = 0.0;
      delta[2] = 0.0;
      deltb[0] = 0.0;
      deltb[1] = 0.0;
      deltb[2] = 0.0;

      // theta1 torque
      if (theta1 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * df4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

        t1dir[0] = d_nx_xtrct(a,1) * d_nx_xtrct(b,2) - d_nx_xtrct(a,2) * d_nx_xtrct(b,1);
        t1dir[1] = d_nx_xtrct(a,2) * d_nx_xtrct(b,0) - d_nx_xtrct(a,0) * d_nx_xtrct(b,2);
        t1dir[2] = d_nx_xtrct(a,0) * d_nx_xtrct(b,1) - d_nx_xtrct(a,1) * d_nx_xtrct(b,0);
        delta[0] += t1dir[0] * tpair;
        delta[1] += t1dir[1] * tpair;
        delta[2] += t1dir[2] * tpair;
        deltb[0] += t1dir[0] * tpair;
        deltb[1] += t1dir[1] * tpair;
        deltb[2] += t1dir[2] * tpair;
      }
      //theta2 torque
      if (theta2 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

        t2dir[0] = d_nx_xtrct(a,1) * delr_hb_norm[2] - d_nx_xtrct(a,2) * delr_hb_norm[1];
        t2dir[1] = d_nx_xtrct(a,2) * delr_hb_norm[0] - d_nx_xtrct(a,0) * delr_hb_norm[2];
        t2dir[2] = d_nx_xtrct(a,0) * delr_hb_norm[1] - d_nx_xtrct(a,1) * delr_hb_norm[0];
        delta[0] += t2dir[0] * tpair;
        delta[1] += t2dir[1] * tpair;
        delta[2] += t2dir[2] * tpair;
      }
      //theta3 torque
      if (theta3 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

        t3dir[0] = d_nx_xtrct(b,1) * delr_hb_norm[2] - d_nx_xtrct(b,2) * delr_hb_norm[1];
        t3dir[1] = d_nx_xtrct(b,2) * delr_hb_norm[0] - d_nx_xtrct(b,0) * delr_hb_norm[2];
        t3dir[2] = d_nx_xtrct(b,0) * delr_hb_norm[1] - d_nx_xtrct(b,1) * delr_hb_norm[0];
        deltb[0] += t3dir[0] * tpair;
        deltb[1] += t3dir[1] * tpair;
        deltb[2] += t3dir[2] * tpair;
      }
      //theta4 torque
      if (theta4 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * df4t4 * f4t7 * f4t8 * factor_lj;

        t4dir[0] = d_nz_xtrct(b,1) * d_nz_xtrct(a,2) - d_nz_xtrct(b,2) * d_nz_xtrct(a,1);
        t4dir[1] = d_nz_xtrct(b,2) * d_nz_xtrct(a,0) - d_nz_xtrct(b,0) * d_nz_xtrct(a,2);
        t4dir[2] = d_nz_xtrct(b,0) * d_nz_xtrct(a,1) - d_nz_xtrct(b,1) * d_nz_xtrct(a,0);
        delta[0] += t4dir[0] * tpair;
        delta[1] += t4dir[1] * tpair;
        delta[2] += t4dir[2] * tpair;
        deltb[0] += t4dir[0] * tpair;
        deltb[1] += t4dir[1] * tpair;
        deltb[2] += t4dir[2] * tpair;
      }
      //theta7 torque
      if (theta7 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * factor_lj;

        t7dir[0] = d_nz_xtrct(a,1) * delr_hb_norm[2] - d_nz_xtrct(a,2) * delr_hb_norm[1];
        t7dir[1] = d_nz_xtrct(a,2) * delr_hb_norm[0] - d_nz_xtrct(a,0) * delr_hb_norm[2];
        t7dir[2] = d_nz_xtrct(a,0) * delr_hb_norm[1] - d_nz_xtrct(a,1) * delr_hb_norm[0];
        delta[0] += t7dir[0] * tpair;
        delta[1] += t7dir[1] * tpair;
        delta[2] += t7dir[2] * tpair;
      }
      //theta8 torque
      if (theta8 != static_cast<KK_FLOAT>(0.0)) {

        tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * factor_lj;

        t8dir[0] = d_nz_xtrct(b,1) * delr_hb_norm[2] - d_nz_xtrct(b,2) * delr_hb_norm[1];
        t8dir[1] = d_nz_xtrct(b,2) * delr_hb_norm[0] - d_nz_xtrct(b,0) * delr_hb_norm[2];
        t8dir[2] = d_nz_xtrct(b,0) * delr_hb_norm[1] - d_nz_xtrct(b,1) * delr_hb_norm[0];
        deltb[0] += t8dir[0] * tpair;
        deltb[1] += t8dir[1] * tpair;
        deltb[2] += t8dir[2] * tpair;
      }

      // increment torques

      a_torque(a,0) += delta[0];
      a_torque(a,1) += delta[1];
      a_torque(a,2) += delta[2];

      if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
        a_torque(b,0) -= deltb[0];
        a_torque(b,1) -= deltb[1];
        a_torque(b,2) -= deltb[2];
      }
    // end of early rejection criterion
    } // evdwl
    } // f4t7
    } // f4t4
    } // f4t3
    } // f4t2
    } // f4t1
    } // f1
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdnaHbondCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ia,ev);
}

/* ----------------------------------------------------------------------
   ComputeGPUPair Functor(s) and staged hbond helpers for lower
   live register pressure in GPU kernels.
-------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_radial_terms(const int &atype, const int &btype,
  const KK_FLOAT &r_hb, KK_FLOAT &f1, KK_FLOAT &df1) const
{
  const KK_FLOAT p_epsilon_hb = d_epsilon_hb(atype,btype);
  const KK_FLOAT p_a_hb = d_a_hb(atype,btype);
  const KK_FLOAT p_cut_hb_0 = d_cut_hb_0(atype,btype);
  const KK_FLOAT p_cut_hb_lc = d_cut_hb_lc(atype,btype);
  const KK_FLOAT p_cut_hb_hc = d_cut_hb_hc(atype,btype);
  const KK_FLOAT p_cut_hb_lo = d_cut_hb_lo(atype,btype);
  const KK_FLOAT p_cut_hb_hi = d_cut_hb_hi(atype,btype);
  const KK_FLOAT p_b_hb_lo = d_b_hb_lo(atype,btype);
  const KK_FLOAT p_b_hb_hi = d_b_hb_hi(atype,btype);
  const KK_FLOAT p_shift_hb = d_shift_hb(atype,btype);

    f1 = F1_KK(r_hb, p_epsilon_hb, p_a_hb, p_cut_hb_0,
      p_cut_hb_lc, p_cut_hb_hc, p_cut_hb_lo, p_cut_hb_hi,
      p_b_hb_lo, p_b_hb_hi, p_shift_hb);
    if (f1 == static_cast<KK_FLOAT>(0.0)) return false;

  df1 = DF1_KK(r_hb, p_epsilon_hb, p_a_hb, p_cut_hb_0,
      p_cut_hb_lc, p_cut_hb_hc, p_cut_hb_lo, p_cut_hb_hi,
      p_b_hb_lo, p_b_hb_hi);
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta1_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  KK_FLOAT &theta1, KK_FLOAT &f4t1, KK_FLOAT &df4t1) const
{
  const KK_FLOAT p_a_hb1 = d_a_hb1(atype,btype);
  const KK_FLOAT p_theta_hb1_0 = d_theta_hb1_0(atype,btype);
  const KK_FLOAT p_dtheta_hb1_ast = d_dtheta_hb1_ast(atype,btype);
  const KK_FLOAT p_b_hb1 = d_b_hb1(atype,btype);
  const KK_FLOAT p_dtheta_hb1_c = d_dtheta_hb1_c(atype,btype);

  KK_FLOAT cost1 = -Kokkos::fma(a_nx[2], b_nx[2], Kokkos::fma(a_nx[1], b_nx[1], a_nx[0] * b_nx[0]));
  if (cost1 > 1.0) cost1 = 1.0;
  if (cost1 < -1.0) cost1 = -1.0;
  theta1 = acos(cost1);

  f4t1 = F4_KK(theta1, p_a_hb1, p_theta_hb1_0, p_dtheta_hb1_ast, p_b_hb1, p_dtheta_hb1_c);
  if (f4t1 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin1_sq = Kokkos::fma(-cost1, cost1, static_cast<KK_FLOAT>(1.0));
  if (sin1_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin1 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin1_sq);
  df4t1 = DF4_KK(theta1, p_a_hb1, p_theta_hb1_0, p_dtheta_hb1_ast, p_b_hb1, p_dtheta_hb1_c) * rsin1;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta2_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta2, KK_FLOAT &cost2, KK_FLOAT &f4t2, KK_FLOAT &df4t2) const
{
  const KK_FLOAT p_a_hb2 = d_a_hb2(atype,btype);
  const KK_FLOAT p_theta_hb2_0 = d_theta_hb2_0(atype,btype);
  const KK_FLOAT p_dtheta_hb2_ast = d_dtheta_hb2_ast(atype,btype);
  const KK_FLOAT p_b_hb2 = d_b_hb2(atype,btype);
  const KK_FLOAT p_dtheta_hb2_c = d_dtheta_hb2_c(atype,btype);

  cost2 = -Kokkos::fma(a_nx[2], delr_hb_norm[2], Kokkos::fma(a_nx[1], delr_hb_norm[1], a_nx[0] * delr_hb_norm[0]));
  if (cost2 > 1.0) cost2 = 1.0;
  if (cost2 < -1.0) cost2 = -1.0;
  theta2 = acos(cost2);

  f4t2 = F4_KK(theta2, p_a_hb2, p_theta_hb2_0, p_dtheta_hb2_ast, p_b_hb2, p_dtheta_hb2_c);
  if (f4t2 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin2_sq = Kokkos::fma(-cost2, cost2, static_cast<KK_FLOAT>(1.0));
  if (sin2_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin2 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin2_sq);
  df4t2 = DF4_KK(theta2, p_a_hb2, p_theta_hb2_0, p_dtheta_hb2_ast, p_b_hb2, p_dtheta_hb2_c) * rsin2;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta3_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta3, KK_FLOAT &cost3, KK_FLOAT &f4t3, KK_FLOAT &df4t3) const
{
  const KK_FLOAT p_a_hb3 = d_a_hb3(atype,btype);
  const KK_FLOAT p_theta_hb3_0 = d_theta_hb3_0(atype,btype);
  const KK_FLOAT p_dtheta_hb3_ast = d_dtheta_hb3_ast(atype,btype);
  const KK_FLOAT p_b_hb3 = d_b_hb3(atype,btype);
  const KK_FLOAT p_dtheta_hb3_c = d_dtheta_hb3_c(atype,btype);

  cost3 = Kokkos::fma(b_nx[2], delr_hb_norm[2], Kokkos::fma(b_nx[1], delr_hb_norm[1], b_nx[0] * delr_hb_norm[0]));
  if (cost3 > 1.0) cost3 = 1.0;
  if (cost3 < -1.0) cost3 = -1.0;
  theta3 = acos(cost3);

  f4t3 = F4_KK(theta3, p_a_hb3, p_theta_hb3_0, p_dtheta_hb3_ast, p_b_hb3, p_dtheta_hb3_c);
  if (f4t3 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin3_sq = Kokkos::fma(-cost3, cost3, static_cast<KK_FLOAT>(1.0));
  if (sin3_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin3 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin3_sq);
  df4t3 = DF4_KK(theta3, p_a_hb3, p_theta_hb3_0, p_dtheta_hb3_ast, p_b_hb3, p_dtheta_hb3_c) * rsin3;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta4_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  KK_FLOAT &theta4, KK_FLOAT &f4t4, KK_FLOAT &df4t4) const
{
  const KK_FLOAT p_a_hb4 = d_a_hb4(atype,btype);
  const KK_FLOAT p_theta_hb4_0 = d_theta_hb4_0(atype,btype);
  const KK_FLOAT p_dtheta_hb4_ast = d_dtheta_hb4_ast(atype,btype);
  const KK_FLOAT p_b_hb4 = d_b_hb4(atype,btype);
  const KK_FLOAT p_dtheta_hb4_c = d_dtheta_hb4_c(atype,btype);

  KK_FLOAT cost4 = Kokkos::fma(a_nz[2], b_nz[2], Kokkos::fma(a_nz[1], b_nz[1], a_nz[0] * b_nz[0]));
  if (cost4 > 1.0) cost4 = 1.0;
  if (cost4 < -1.0) cost4 = -1.0;
  theta4 = acos(cost4);

  f4t4 = F4_KK(theta4, p_a_hb4, p_theta_hb4_0, p_dtheta_hb4_ast, p_b_hb4, p_dtheta_hb4_c);
  if (f4t4 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin4_sq = Kokkos::fma(-cost4, cost4, static_cast<KK_FLOAT>(1.0));
  if (sin4_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin4 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin4_sq);
  df4t4 = DF4_KK(theta4, p_a_hb4, p_theta_hb4_0, p_dtheta_hb4_ast, p_b_hb4, p_dtheta_hb4_c) * rsin4;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta7_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta7, KK_FLOAT &cost7, KK_FLOAT &f4t7, KK_FLOAT &df4t7) const
{
  const KK_FLOAT p_a_hb7 = d_a_hb7(atype,btype);
  const KK_FLOAT p_theta_hb7_0 = d_theta_hb7_0(atype,btype);
  const KK_FLOAT p_dtheta_hb7_ast = d_dtheta_hb7_ast(atype,btype);
  const KK_FLOAT p_b_hb7 = d_b_hb7(atype,btype);
  const KK_FLOAT p_dtheta_hb7_c = d_dtheta_hb7_c(atype,btype);

  cost7 = -Kokkos::fma(a_nz[2], delr_hb_norm[2], Kokkos::fma(a_nz[1], delr_hb_norm[1], a_nz[0] * delr_hb_norm[0]));
  if (cost7 > 1.0) cost7 = 1.0;
  if (cost7 < -1.0) cost7 = -1.0;
  theta7 = acos(cost7);

  f4t7 = F4_KK(theta7, p_a_hb7, p_theta_hb7_0, p_dtheta_hb7_ast, p_b_hb7, p_dtheta_hb7_c);
  if (f4t7 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin7_sq = Kokkos::fma(-cost7, cost7, static_cast<KK_FLOAT>(1.0));
  if (sin7_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin7 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin7_sq);
  df4t7 = DF4_KK(theta7, p_a_hb7, p_theta_hb7_0, p_dtheta_hb7_ast, p_b_hb7, p_dtheta_hb7_c) * rsin7;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta8_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta8, KK_FLOAT &cost8, KK_FLOAT &f4t8, KK_FLOAT &df4t8) const
{
  const KK_FLOAT p_a_hb8 = d_a_hb8(atype,btype);
  const KK_FLOAT p_theta_hb8_0 = d_theta_hb8_0(atype,btype);
  const KK_FLOAT p_dtheta_hb8_ast = d_dtheta_hb8_ast(atype,btype);
  const KK_FLOAT p_b_hb8 = d_b_hb8(atype,btype);
  const KK_FLOAT p_dtheta_hb8_c = d_dtheta_hb8_c(atype,btype);

  cost8 = Kokkos::fma(b_nz[2], delr_hb_norm[2], Kokkos::fma(b_nz[1], delr_hb_norm[1], b_nz[0] * delr_hb_norm[0]));
  if (cost8 > 1.0) cost8 = 1.0;
  if (cost8 < -1.0) cost8 = -1.0;
  theta8 = acos(cost8);

  f4t8 = F4_KK(theta8, p_a_hb8, p_theta_hb8_0, p_dtheta_hb8_ast, p_b_hb8, p_dtheta_hb8_c);
  if (f4t8 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin8_sq = Kokkos::fma(-cost8, cost8, static_cast<KK_FLOAT>(1.0));
  if (sin8_sq <= static_cast<KK_FLOAT>(0.0)) return false;
  const KK_FLOAT rsin8 = static_cast<KK_FLOAT>(1.0) / sqrtf(sin8_sq);
  df4t8 = DF4_KK(theta8, p_a_hb8, p_theta_hb8_0, p_dtheta_hb8_ast, p_b_hb8, p_dtheta_hb8_c) * rsin8;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::hbond_force_contrib(const KK_FLOAT &f1,
  const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
  const KK_FLOAT &f4t4, const KK_FLOAT &f4t7, const KK_FLOAT &f4t8,
  const KK_FLOAT &df1, const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
  const KK_FLOAT &df4t7, const KK_FLOAT &df4t8,
  const KK_FLOAT &rinv_hb, const KK_FLOAT &factor_lj,
  const KK_FLOAT &theta2, const KK_FLOAT &theta3, const KK_FLOAT &theta7, const KK_FLOAT &theta8,
  const KK_FLOAT &cost2, const KK_FLOAT &cost3, const KK_FLOAT &cost7, const KK_FLOAT &cost8,
  const KK_FLOAT (&delr_hb)[3], const KK_FLOAT (&delr_hb_norm)[3],
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&ra_chb)[3], const KK_FLOAT (&rb_chb)[3],
  KK_ACC_FLOAT (&delf)[3], KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  KK_FLOAT finc = -df1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

  delf[0] = Kokkos::fma(delr_hb[0], finc, delf[0]);
  delf[1] = Kokkos::fma(delr_hb[1], finc, delf[1]);
  delf[2] = Kokkos::fma(delr_hb[2], finc, delf[2]);

  if (theta2 != static_cast<KK_FLOAT>(0.0)) {
    finc = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
    const KK_FLOAT t2f0 = Kokkos::fma(delr_hb_norm[0], cost2, a_nx[0]);
    const KK_FLOAT t2f1 = Kokkos::fma(delr_hb_norm[1], cost2, a_nx[1]);
    const KK_FLOAT t2f2 = Kokkos::fma(delr_hb_norm[2], cost2, a_nx[2]);
    delf[0] = Kokkos::fma(t2f0, finc, delf[0]);
    delf[1] = Kokkos::fma(t2f1, finc, delf[1]);
    delf[2] = Kokkos::fma(t2f2, finc, delf[2]);
  }

  if (theta3 != static_cast<KK_FLOAT>(0.0)) {
    finc = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj;
    const KK_FLOAT t3f0 = Kokkos::fma(delr_hb_norm[0], cost3, -b_nx[0]);
    const KK_FLOAT t3f1 = Kokkos::fma(delr_hb_norm[1], cost3, -b_nx[1]);
    const KK_FLOAT t3f2 = Kokkos::fma(delr_hb_norm[2], cost3, -b_nx[2]);
    delf[0] = Kokkos::fma(t3f0, finc, delf[0]);
    delf[1] = Kokkos::fma(t3f1, finc, delf[1]);
    delf[2] = Kokkos::fma(t3f2, finc, delf[2]);
  }

  if (theta7 != static_cast<KK_FLOAT>(0.0)) {
    finc = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * rinv_hb * factor_lj;
    const KK_FLOAT t7f0 = Kokkos::fma(delr_hb_norm[0], cost7, a_nz[0]);
    const KK_FLOAT t7f1 = Kokkos::fma(delr_hb_norm[1], cost7, a_nz[1]);
    const KK_FLOAT t7f2 = Kokkos::fma(delr_hb_norm[2], cost7, a_nz[2]);
    delf[0] = Kokkos::fma(t7f0, finc, delf[0]);
    delf[1] = Kokkos::fma(t7f1, finc, delf[1]);
    delf[2] = Kokkos::fma(t7f2, finc, delf[2]);
  }

  if (theta8 != static_cast<KK_FLOAT>(0.0)) {
    finc = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * rinv_hb * factor_lj;
    const KK_FLOAT t8f0 = Kokkos::fma(delr_hb_norm[0], cost8, -b_nz[0]);
    const KK_FLOAT t8f1 = Kokkos::fma(delr_hb_norm[1], cost8, -b_nz[1]);
    const KK_FLOAT t8f2 = Kokkos::fma(delr_hb_norm[2], cost8, -b_nz[2]);
    delf[0] = Kokkos::fma(t8f0, finc, delf[0]);
    delf[1] = Kokkos::fma(t8f1, finc, delf[1]);
    delf[2] = Kokkos::fma(t8f2, finc, delf[2]);
  }

  delta[0] = Kokkos::fma(ra_chb[1], delf[2], -ra_chb[2] * delf[1]);
  delta[1] = Kokkos::fma(ra_chb[2], delf[0], -ra_chb[0] * delf[2]);
  delta[2] = Kokkos::fma(ra_chb[0], delf[1], -ra_chb[1] * delf[0]);

  deltb[0] = Kokkos::fma(rb_chb[1], delf[2], -rb_chb[2] * delf[1]);
  deltb[1] = Kokkos::fma(rb_chb[2], delf[0], -rb_chb[0] * delf[2]);
  deltb[2] = Kokkos::fma(rb_chb[0], delf[1], -rb_chb[1] * delf[0]);
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::hbond_torque_contrib(const KK_FLOAT &f1,
  const KK_FLOAT &f4t1, const KK_FLOAT &f4t2, const KK_FLOAT &f4t3,
  const KK_FLOAT &f4t4, const KK_FLOAT &f4t7, const KK_FLOAT &f4t8,
  const KK_FLOAT &df4t1, const KK_FLOAT &df4t2, const KK_FLOAT &df4t3,
  const KK_FLOAT &df4t4, const KK_FLOAT &df4t7, const KK_FLOAT &df4t8,
  const KK_FLOAT &factor_lj,
  const KK_FLOAT &theta1, const KK_FLOAT &theta2, const KK_FLOAT &theta3,
  const KK_FLOAT &theta4, const KK_FLOAT &theta7, const KK_FLOAT &theta8,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&b_nx)[3],
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  const KK_FLOAT (&delr_hb_norm)[3],
  KK_ACC_FLOAT (&delta)[3], KK_ACC_FLOAT (&deltb)[3]) const
{
  delta[0] = 0.0;
  delta[1] = 0.0;
  delta[2] = 0.0;
  deltb[0] = 0.0;
  deltb[1] = 0.0;
  deltb[2] = 0.0;

  KK_FLOAT tpair;

  if (theta1 != static_cast<KK_FLOAT>(0.0)) {
    tpair = -f1 * df4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

    const KK_FLOAT t1dir0 = Kokkos::fma(a_nx[1], b_nx[2], -a_nx[2] * b_nx[1]);
    const KK_FLOAT t1dir1 = Kokkos::fma(a_nx[2], b_nx[0], -a_nx[0] * b_nx[2]);
    const KK_FLOAT t1dir2 = Kokkos::fma(a_nx[0], b_nx[1], -a_nx[1] * b_nx[0]);
    delta[0] = Kokkos::fma(t1dir0, tpair, delta[0]);
    delta[1] = Kokkos::fma(t1dir1, tpair, delta[1]);
    delta[2] = Kokkos::fma(t1dir2, tpair, delta[2]);
    deltb[0] = Kokkos::fma(t1dir0, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t1dir1, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t1dir2, tpair, deltb[2]);
  }

  if (theta2 != static_cast<KK_FLOAT>(0.0)) {
    tpair = -f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

    const KK_FLOAT t2dir0 = Kokkos::fma(a_nx[1], delr_hb_norm[2], -a_nx[2] * delr_hb_norm[1]);
    const KK_FLOAT t2dir1 = Kokkos::fma(a_nx[2], delr_hb_norm[0], -a_nx[0] * delr_hb_norm[2]);
    const KK_FLOAT t2dir2 = Kokkos::fma(a_nx[0], delr_hb_norm[1], -a_nx[1] * delr_hb_norm[0]);
    delta[0] = Kokkos::fma(t2dir0, tpair, delta[0]);
    delta[1] = Kokkos::fma(t2dir1, tpair, delta[1]);
    delta[2] = Kokkos::fma(t2dir2, tpair, delta[2]);
  }

  if (theta3 != static_cast<KK_FLOAT>(0.0)) {
    tpair = -f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * factor_lj;

    const KK_FLOAT t3dir0 = Kokkos::fma(b_nx[1], delr_hb_norm[2], -b_nx[2] * delr_hb_norm[1]);
    const KK_FLOAT t3dir1 = Kokkos::fma(b_nx[2], delr_hb_norm[0], -b_nx[0] * delr_hb_norm[2]);
    const KK_FLOAT t3dir2 = Kokkos::fma(b_nx[0], delr_hb_norm[1], -b_nx[1] * delr_hb_norm[0]);
    deltb[0] = Kokkos::fma(t3dir0, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t3dir1, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t3dir2, tpair, deltb[2]);
  }

  if (theta4) {
    tpair = -f1 * f4t1 * f4t2 * f4t3 * df4t4 * f4t7 * f4t8 * factor_lj;

    const KK_FLOAT t4dir0 = Kokkos::fma(b_nz[1], a_nz[2], -b_nz[2] * a_nz[1]);
    const KK_FLOAT t4dir1 = Kokkos::fma(b_nz[2], a_nz[0], -b_nz[0] * a_nz[2]);
    const KK_FLOAT t4dir2 = Kokkos::fma(b_nz[0], a_nz[1], -b_nz[1] * a_nz[0]);
    delta[0] = Kokkos::fma(t4dir0, tpair, delta[0]);
    delta[1] = Kokkos::fma(t4dir1, tpair, delta[1]);
    delta[2] = Kokkos::fma(t4dir2, tpair, delta[2]);
    deltb[0] = Kokkos::fma(t4dir0, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t4dir1, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t4dir2, tpair, deltb[2]);
  }

  if (theta7 != static_cast<KK_FLOAT>(0.0)) {
    tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * factor_lj;

    const KK_FLOAT t7dir0 = Kokkos::fma(a_nz[1], delr_hb_norm[2], -a_nz[2] * delr_hb_norm[1]);
    const KK_FLOAT t7dir1 = Kokkos::fma(a_nz[2], delr_hb_norm[0], -a_nz[0] * delr_hb_norm[2]);
    const KK_FLOAT t7dir2 = Kokkos::fma(a_nz[0], delr_hb_norm[1], -a_nz[1] * delr_hb_norm[0]);
    delta[0] = Kokkos::fma(t7dir0, tpair, delta[0]);
    delta[1] = Kokkos::fma(t7dir1, tpair, delta[1]);
    delta[2] = Kokkos::fma(t7dir2, tpair, delta[2]);
  }

  if (theta8 != static_cast<KK_FLOAT>(0.0)) {
    tpair = -f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * factor_lj;

    const KK_FLOAT t8dir0 = Kokkos::fma(b_nz[1], delr_hb_norm[2], -b_nz[2] * delr_hb_norm[1]);
    const KK_FLOAT t8dir1 = Kokkos::fma(b_nz[2], delr_hb_norm[0], -b_nz[0] * delr_hb_norm[2]);
    const KK_FLOAT t8dir2 = Kokkos::fma(b_nz[0], delr_hb_norm[1], -b_nz[1] * delr_hb_norm[0]);
    deltb[0] = Kokkos::fma(t8dir0, tpair, deltb[0]);
    deltb[1] = Kokkos::fma(t8dir1, tpair, deltb[1]);
    deltb[2] = Kokkos::fma(t8dir2, tpair, deltb[2]);
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair, EV_FLOAT &ev) const
{
  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // Direct packed pair lookup: high 32 bits = a, low 32 bits = b.
  const uint64_t pair = d_pairs_screened(ipair);
  // "pair >> 32" shifts the pair to the right by 32 bits, so the upper 32 bits
  // becomes the lower 32 bits to recover the atom-a index.
  const int a = static_cast<int>(pair >> 32);
  const int atype = type(a);
  // "pair & 0xffffffffu" keeps only the lower 32 bits to recover the atom-b index.
  int b = static_cast<int>(pair & 0xffffffffu);
  const KK_FLOAT factor_lj = special_lj[sbmask(b)];
  if (factor_lj == static_cast<KK_FLOAT>(0.0)) return;
  b &= NEIGHMASK;
  const int btype = type(b);

  if (unique_basepair_enabled) {
    const int idca = d_idc(a);
    const int idcb = d_idc(b);
    if (idca != tag(b) && idcb != tag(a) && idca > 0 && idcb > 0) return;
  }

  KK_FLOAT a_nx[3], a_nz[3], b_nx[3], b_nz[3];
  KK_FLOAT ra_chb[3], rb_chb[3];
  KK_FLOAT delr_hb[3], delr_hb_norm[3];
  KK_FLOAT rsq_hb, r_hb, rinv_hb;
  KK_FLOAT evdwl;

  a_nx[0] = d_nx_xtrct(a,0);
  a_nx[1] = d_nx_xtrct(a,1);
  a_nx[2] = d_nx_xtrct(a,2);
  a_nz[0] = d_nz_xtrct(a,0);
  a_nz[1] = d_nz_xtrct(a,1);
  a_nz[2] = d_nz_xtrct(a,2);

  b_nx[0] = d_nx_xtrct(b,0);
  b_nx[1] = d_nx_xtrct(b,1);
  b_nx[2] = d_nx_xtrct(b,2);
  b_nz[0] = d_nz_xtrct(b,0);
  b_nz[1] = d_nz_xtrct(b,1);
  b_nz[2] = d_nz_xtrct(b,2);

  // vector COM-hbond site a and b
  if (OXDNAFLAG == OXDNA) {
    // Used for oxdna1 and oxdna2, which have the same hbond site offset.
    constexpr KK_FLOAT dx_cbs_oxdna1=+0.4;
    ra_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
    rb_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(b,0);
    rb_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(b,1);
    rb_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(b,2);
  } else if (OXDNAFLAG == OXDNA3) {
    constexpr KK_FLOAT dx_cbs_pur_oxdna3 = +0.43;
    constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = +0.37;
    int nucl_acid = (atype%4);
    if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
      ra_chb[0] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,0);
      ra_chb[1] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,1);
      ra_chb[2] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,2);
    } else {  // purine (A or G)
      ra_chb[0] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,0);
      ra_chb[1] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,1);
      ra_chb[2] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,2);
    }
    nucl_acid = (btype%4);
    if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
      rb_chb[0] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,0);
      rb_chb[1] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,1);
      rb_chb[2] = dx_cbs_pyr_oxdna3*d_nx_xtrct(b,2);
    } else {  // purine (A or G)
      rb_chb[0] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,0);
      rb_chb[1] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,1);
      rb_chb[2] = dx_cbs_pur_oxdna3*d_nx_xtrct(b,2);
    }
  }

  delr_hb[0] = x(a,0) + ra_chb[0] - x(b,0) - rb_chb[0];
  delr_hb[1] = x(a,1) + ra_chb[1] - x(b,1) - rb_chb[1];
  delr_hb[2] = x(a,2) + ra_chb[2] - x(b,2) - rb_chb[2];

  rsq_hb = Kokkos::fma(delr_hb[2], delr_hb[2],
      Kokkos::fma(delr_hb[1], delr_hb[1], delr_hb[0] * delr_hb[0]));
  if (rsq_hb <= static_cast<KK_FLOAT>(0.0)) return;
  rinv_hb = static_cast<KK_FLOAT>(1.0) / sqrtf(rsq_hb);
  r_hb = rsq_hb * rinv_hb;

  delr_hb_norm[0] = delr_hb[0] * rinv_hb;
  delr_hb_norm[1] = delr_hb[1] * rinv_hb;
  delr_hb_norm[2] = delr_hb[2] * rinv_hb;

  KK_FLOAT f1, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8;
  KK_FLOAT df1, df4t1, df4t2, df4t3, df4t4, df4t7, df4t8;
  KK_FLOAT theta1, theta2, theta3, theta4, theta7, theta8;
  KK_FLOAT cost2, cost3, cost7, cost8;

  if (!hbond_radial_terms(atype, btype, r_hb, f1, df1)) return;
  if (!hbond_theta1_terms(atype, btype, a_nx, b_nx, theta1, f4t1, df4t1)) return;
  if (!hbond_theta2_terms(atype, btype, a_nx, delr_hb_norm, theta2, cost2, f4t2, df4t2)) return;
  if (!hbond_theta3_terms(atype, btype, b_nx, delr_hb_norm, theta3, cost3, f4t3, df4t3)) return;
  if (!hbond_theta4_terms(atype, btype, a_nz, b_nz, theta4, f4t4, df4t4)) return;
  if (!hbond_theta7_terms(atype, btype, a_nz, delr_hb_norm, theta7, cost7, f4t7, df4t7)) return;
  if (!hbond_theta8_terms(atype, btype, b_nz, delr_hb_norm, theta8, cost8, f4t8, df4t8)) return;

  evdwl = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
  if (evdwl == static_cast<KK_FLOAT>(0.0)) return;

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];
  delf[0] = 0.0;
  delf[1] = 0.0;
  delf[2] = 0.0;
  delta[0] = 0.0;
  delta[1] = 0.0;
  delta[2] = 0.0;
  deltb[0] = 0.0;
  deltb[1] = 0.0;
  deltb[2] = 0.0;

  hbond_force_contrib(
    f1, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8,
    df1, df4t2, df4t3, df4t7, df4t8,
    rinv_hb, factor_lj,
    theta2, theta3, theta7, theta8,
    cost2, cost3, cost7, cost8,
    delr_hb, delr_hb_norm,
    a_nx, b_nx, a_nz, b_nz,
    ra_chb, rb_chb,
    delf, delta, deltb);

  a_f(a,0) += delf[0];
  a_f(a,1) += delf[1];
  a_f(a,2) += delf[2];
  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_f(b,0) -= delf[0];
    a_f(b,1) -= delf[1];
    a_f(b,2) -= delf[2];
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }

  if (EVFLAG) {
    ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;

    if (vflag_either || eflag_atom) {
      this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl,\
      delf[0],delf[1],delf[2],x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
    }
  }

  hbond_torque_contrib(
    f1, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8,
    df4t1, df4t2, df4t3, df4t4, df4t7, df4t8,
    factor_lj,
    theta1, theta2, theta3, theta4, theta7, theta8,
    a_nx, b_nx, a_nz, b_nz, delr_hb_norm,
    delta, deltb);

  a_torque(a,0) += delta[0];
  a_torque(a,1) += delta[1];
  a_torque(a,2) += delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_torque(b,0) -= deltb[0];
    a_torque(b,1) -= deltb[1];
    a_torque(b,2) -= deltb[2];
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ipair,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::allocate()
{
  PairOxdnaHbond::allocate();

  int n = atom->ntypes;

  memoryKK->create_kokkos(k_epsilon_hb,n+1,n+1,"PairOxdnaHbond:epsilon_hb");
  memoryKK->create_kokkos(k_a_hb,n+1,n+1,"PairOxdnaHbond:a_hb");
  memoryKK->create_kokkos(k_cut_hb_0,n+1,n+1,"PairOxdnaHbond:cut_hb_0");
  memoryKK->create_kokkos(k_cut_hb_c,n+1,n+1,"PairOxdnaHbond:cut_hb_c");
  memoryKK->create_kokkos(k_cut_hb_lo,n+1,n+1,"PairOxdnaHbond:cut_hb_lo");
  memoryKK->create_kokkos(k_cut_hb_hi,n+1,n+1,"PairOxdnaHbond:cut_hb_hi");
  memoryKK->create_kokkos(k_cut_hb_lc,n+1,n+1,"PairOxdnaHbond:cut_hb_lc");
  memoryKK->create_kokkos(k_cut_hb_hc,n+1,n+1,"PairOxdnaHbond:cut_hb_hc");
  memoryKK->create_kokkos(k_b_hb_lo,n+1,n+1,"PairOxdnaHbond:b_hb_lo");
  memoryKK->create_kokkos(k_b_hb_hi,n+1,n+1,"PairOxdnaHbond:b_hb_hi");
  memoryKK->create_kokkos(k_shift_hb,n+1,n+1,"PairOxdnaHbond:shift_hb");
  memoryKK->create_kokkos(k_cutsq_hb_hc,n+1,n+1,"PairOxdnaHbond:cutsq_hb_hc");

  memoryKK->create_kokkos(k_a_hb1,n+1,n+1,"PairOxdnaHbond:a_hb1");
  memoryKK->create_kokkos(k_theta_hb1_0,n+1,n+1,"PairOxdnaHbond:theta_hb1_0");
  memoryKK->create_kokkos(k_dtheta_hb1_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb1_ast");
  memoryKK->create_kokkos(k_b_hb1,n+1,n+1,"PairOxdnaHbond:b_hb1");
  memoryKK->create_kokkos(k_dtheta_hb1_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb1_c");

  memoryKK->create_kokkos(k_a_hb2,n+1,n+1,"PairOxdnaHbond:a_hb2");
  memoryKK->create_kokkos(k_theta_hb2_0,n+1,n+1,"PairOxdnaHbond:theta_hb2_0");
  memoryKK->create_kokkos(k_dtheta_hb2_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb2_ast");
  memoryKK->create_kokkos(k_b_hb2,n+1,n+1,"PairOxdnaHbond:b_hb2");
  memoryKK->create_kokkos(k_dtheta_hb2_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb2_c");

  memoryKK->create_kokkos(k_a_hb3,n+1,n+1,"PairOxdnaHbond:a_hb3");
  memoryKK->create_kokkos(k_theta_hb3_0,n+1,n+1,"PairOxdnaHbond:theta_hb3_0");
  memoryKK->create_kokkos(k_dtheta_hb3_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb3_ast");
  memoryKK->create_kokkos(k_b_hb3,n+1,n+1,"PairOxdnaHbond:b_hb3");
  memoryKK->create_kokkos(k_dtheta_hb3_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb3_c");

  memoryKK->create_kokkos(k_a_hb4,n+1,n+1,"PairOxdnaHbond:a_hb4");
  memoryKK->create_kokkos(k_theta_hb4_0,n+1,n+1,"PairOxdnaHbond:theta_hb4_0");
  memoryKK->create_kokkos(k_dtheta_hb4_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb4_ast");
  memoryKK->create_kokkos(k_b_hb4,n+1,n+1,"PairOxdnaHbond:b_hb4");
  memoryKK->create_kokkos(k_dtheta_hb4_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb4_c");

  memoryKK->create_kokkos(k_a_hb7,n+1,n+1,"PairOxdnaHbond:a_hb7");
  memoryKK->create_kokkos(k_theta_hb7_0,n+1,n+1,"PairOxdnaHbond:theta_hb7_0");
  memoryKK->create_kokkos(k_dtheta_hb7_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb7_ast");
  memoryKK->create_kokkos(k_b_hb7,n+1,n+1,"PairOxdnaHbond:b_hb7");
  memoryKK->create_kokkos(k_dtheta_hb7_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb7_c");

  memoryKK->create_kokkos(k_a_hb8,n+1,n+1,"PairOxdnaHbond:a_hb8");
  memoryKK->create_kokkos(k_theta_hb8_0,n+1,n+1,"PairOxdnaHbond:theta_hb8_0");
  memoryKK->create_kokkos(k_dtheta_hb8_ast,n+1,n+1,"PairOxdnaHbond:dtheta_hb8_ast");
  memoryKK->create_kokkos(k_b_hb8,n+1,n+1,"PairOxdnaHbond:b_hb8");
  memoryKK->create_kokkos(k_dtheta_hb8_c,n+1,n+1,"PairOxdnaHbond:dtheta_hb8_c");

  d_epsilon_hb = k_epsilon_hb.template view<DeviceType>();
  d_a_hb = k_a_hb.template view<DeviceType>();
  d_cut_hb_0 = k_cut_hb_0.template view<DeviceType>();
  d_cut_hb_c = k_cut_hb_c.template view<DeviceType>();
  d_cut_hb_lo = k_cut_hb_lo.template view<DeviceType>();
  d_cut_hb_hi = k_cut_hb_hi.template view<DeviceType>();
  d_cut_hb_lc = k_cut_hb_lc.template view<DeviceType>();
  d_cut_hb_hc = k_cut_hb_hc.template view<DeviceType>();
  d_b_hb_lo = k_b_hb_lo.template view<DeviceType>();
  d_b_hb_hi = k_b_hb_hi.template view<DeviceType>();
  d_shift_hb = k_shift_hb.template view<DeviceType>();
  d_cutsq_hb_hc = k_cutsq_hb_hc.template view<DeviceType>();

  d_a_hb1 = k_a_hb1.template view<DeviceType>();
  d_theta_hb1_0 = k_theta_hb1_0.template view<DeviceType>();
  d_dtheta_hb1_ast = k_dtheta_hb1_ast.template view<DeviceType>();
  d_b_hb1 = k_b_hb1.template view<DeviceType>();
  d_dtheta_hb1_c = k_dtheta_hb1_c.template view<DeviceType>();

  d_a_hb2 = k_a_hb2.template view<DeviceType>();
  d_theta_hb2_0 = k_theta_hb2_0.template view<DeviceType>();
  d_dtheta_hb2_ast = k_dtheta_hb2_ast.template view<DeviceType>();
  d_b_hb2 = k_b_hb2.template view<DeviceType>();
  d_dtheta_hb2_c = k_dtheta_hb2_c.template view<DeviceType>();

  d_a_hb3 = k_a_hb3.template view<DeviceType>();
  d_theta_hb3_0 = k_theta_hb3_0.template view<DeviceType>();
  d_dtheta_hb3_ast = k_dtheta_hb3_ast.template view<DeviceType>();
  d_b_hb3 = k_b_hb3.template view<DeviceType>();
  d_dtheta_hb3_c = k_dtheta_hb3_c.template view<DeviceType>();

  d_a_hb4 = k_a_hb4.template view<DeviceType>();
  d_theta_hb4_0 = k_theta_hb4_0.template view<DeviceType>();
  d_dtheta_hb4_ast = k_dtheta_hb4_ast.template view<DeviceType>();
  d_b_hb4 = k_b_hb4.template view<DeviceType>();
  d_dtheta_hb4_c = k_dtheta_hb4_c.template view<DeviceType>();

  d_a_hb7 = k_a_hb7.template view<DeviceType>();
  d_theta_hb7_0 = k_theta_hb7_0.template view<DeviceType>();
  d_dtheta_hb7_ast = k_dtheta_hb7_ast.template view<DeviceType>();
  d_b_hb7 = k_b_hb7.template view<DeviceType>();
  d_dtheta_hb7_c = k_dtheta_hb7_c.template view<DeviceType>();

  d_a_hb8 = k_a_hb8.template view<DeviceType>();
  d_theta_hb8_0 = k_theta_hb8_0.template view<DeviceType>();
  d_dtheta_hb8_ast = k_dtheta_hb8_ast.template view<DeviceType>();
  d_b_hb8 = k_b_hb8.template view<DeviceType>();
  d_dtheta_hb8_c = k_dtheta_hb8_c.template view<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::init_style()
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

  unique_basepair_enabled = 0;
  idc = nullptr;
  const int ifix = modify->find_fix("Basepairs");
  if (ifix >= 0) {
    int idx, flag, cols;
    idx = atom->find_custom("idc", flag, cols);
    if (idx >= 0 && flag == 0) {
      idc = atom->ivector[idx];
      unique_basepair_enabled = 1;
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdnaHbondKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdnaHbond::init_one(i,j);

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  k_epsilon_hb.view_host()(i,j) = epsilon_hb[i][j]; k_epsilon_hb.view_host()(j,i) = epsilon_hb[j][i];
  k_a_hb.view_host()(i,j) = a_hb[i][j]; k_a_hb.view_host()(j,i) = a_hb[j][i];
  k_cut_hb_0.view_host()(i,j) = cut_hb_0[i][j]; k_cut_hb_0.view_host()(j,i) = cut_hb_0[j][i];
  k_cut_hb_c.view_host()(i,j) = cut_hb_c[i][j]; k_cut_hb_c.view_host()(j,i) = cut_hb_c[j][i];
  k_cut_hb_lo.view_host()(i,j) = cut_hb_lo[i][j]; k_cut_hb_lo.view_host()(j,i) = cut_hb_lo[j][i];
  k_cut_hb_hi.view_host()(i,j) = cut_hb_hi[i][j]; k_cut_hb_hi.view_host()(j,i) = cut_hb_hi[j][i];
  k_cut_hb_lc.view_host()(i,j) = cut_hb_lc[i][j]; k_cut_hb_lc.view_host()(j,i) = cut_hb_lc[j][i];
  k_cut_hb_hc.view_host()(i,j) = cut_hb_hc[i][j]; k_cut_hb_hc.view_host()(j,i) = cut_hb_hc[j][i];
  k_b_hb_lo.view_host()(i,j) = b_hb_lo[i][j]; k_b_hb_lo.view_host()(j,i) = b_hb_lo[j][i];
  k_b_hb_hi.view_host()(i,j) = b_hb_hi[i][j]; k_b_hb_hi.view_host()(j,i) = b_hb_hi[j][i];
  k_shift_hb.view_host()(i,j) = shift_hb[i][j]; k_shift_hb.view_host()(j,i) = shift_hb[j][i];
  k_cutsq_hb_hc.view_host()(i,j) = cutsq_hb_hc[i][j]; k_cutsq_hb_hc.view_host()(j,i) = cutsq_hb_hc[j][i];

  k_a_hb1.view_host()(i,j) = a_hb1[i][j]; k_a_hb1.view_host()(j,i) = a_hb1[j][i];
  k_theta_hb1_0.view_host()(i,j) = theta_hb1_0[i][j]; k_theta_hb1_0.view_host()(j,i) = theta_hb1_0[j][i];
  k_dtheta_hb1_ast.view_host()(i,j) = dtheta_hb1_ast[i][j]; k_dtheta_hb1_ast.view_host()(j,i) = dtheta_hb1_ast[j][i];
  k_b_hb1.view_host()(i,j) = b_hb1[i][j]; k_b_hb1.view_host()(j,i) = b_hb1[j][i];
  k_dtheta_hb1_c.view_host()(i,j) = dtheta_hb1_c[i][j]; k_dtheta_hb1_c.view_host()(j,i) = dtheta_hb1_c[j][i];

  k_a_hb2.view_host()(i,j) = a_hb2[i][j]; k_a_hb2.view_host()(j,i) = a_hb2[j][i];
  k_theta_hb2_0.view_host()(i,j) = theta_hb2_0[i][j]; k_theta_hb2_0.view_host()(j,i) = theta_hb2_0[j][i];
  k_dtheta_hb2_ast.view_host()(i,j) = dtheta_hb2_ast[i][j]; k_dtheta_hb2_ast.view_host()(j,i) = dtheta_hb2_ast[j][i];
  k_b_hb2.view_host()(i,j) = b_hb2[i][j]; k_b_hb2.view_host()(j,i) = b_hb2[j][i];
  k_dtheta_hb2_c.view_host()(i,j) = dtheta_hb2_c[i][j]; k_dtheta_hb2_c.view_host()(j,i) = dtheta_hb2_c[j][i];

  k_a_hb3.view_host()(i,j) = a_hb3[i][j]; k_a_hb3.view_host()(j,i) = a_hb3[j][i];
  k_theta_hb3_0.view_host()(i,j) = theta_hb3_0[i][j]; k_theta_hb3_0.view_host()(j,i) = theta_hb3_0[j][i];
  k_dtheta_hb3_ast.view_host()(i,j) = dtheta_hb3_ast[i][j]; k_dtheta_hb3_ast.view_host()(j,i) = dtheta_hb3_ast[j][i];
  k_b_hb3.view_host()(i,j) = b_hb3[i][j]; k_b_hb3.view_host()(j,i) = b_hb3[j][i];
  k_dtheta_hb3_c.view_host()(i,j) = dtheta_hb3_c[i][j]; k_dtheta_hb3_c.view_host()(j,i) = dtheta_hb3_c[j][i];

  k_a_hb4.view_host()(i,j) = a_hb4[i][j]; k_a_hb4.view_host()(j,i) = a_hb4[j][i];
  k_theta_hb4_0.view_host()(i,j) = theta_hb4_0[i][j]; k_theta_hb4_0.view_host()(j,i) = theta_hb4_0[j][i];
  k_dtheta_hb4_ast.view_host()(i,j) = dtheta_hb4_ast[i][j]; k_dtheta_hb4_ast.view_host()(j,i) = dtheta_hb4_ast[j][i];
  k_b_hb4.view_host()(i,j) = b_hb4[i][j]; k_b_hb4.view_host()(j,i) = b_hb4[j][i];
  k_dtheta_hb4_c.view_host()(i,j) = dtheta_hb4_c[i][j]; k_dtheta_hb4_c.view_host()(j,i) = dtheta_hb4_c[j][i];

  k_a_hb7.view_host()(i,j) = a_hb7[i][j]; k_a_hb7.view_host()(j,i) = a_hb7[j][i];
  k_theta_hb7_0.view_host()(i,j) = theta_hb7_0[i][j]; k_theta_hb7_0.view_host()(j,i) = theta_hb7_0[j][i];
  k_dtheta_hb7_ast.view_host()(i,j) = dtheta_hb7_ast[i][j]; k_dtheta_hb7_ast.view_host()(j,i) = dtheta_hb7_ast[j][i];
  k_b_hb7.view_host()(i,j) = b_hb7[i][j]; k_b_hb7.view_host()(j,i) = b_hb7[j][i];
  k_dtheta_hb7_c.view_host()(i,j) = dtheta_hb7_c[i][j]; k_dtheta_hb7_c.view_host()(j,i) = dtheta_hb7_c[j][i];

  k_a_hb8.view_host()(i,j) = a_hb8[i][j]; k_a_hb8.view_host()(j,i) = a_hb8[j][i];
  k_theta_hb8_0.view_host()(i,j) = theta_hb8_0[i][j]; k_theta_hb8_0.view_host()(j,i) = theta_hb8_0[j][i];
  k_dtheta_hb8_ast.view_host()(i,j) = dtheta_hb8_ast[i][j]; k_dtheta_hb8_ast.view_host()(j,i) = dtheta_hb8_ast[j][i];
  k_b_hb8.view_host()(i,j) = b_hb8[i][j]; k_b_hb8.view_host()(j,i) = b_hb8[j][i];
  k_dtheta_hb8_c.view_host()(i,j) = dtheta_hb8_c[i][j]; k_dtheta_hb8_c.view_host()(j,i) = dtheta_hb8_c[j][i];

  k_epsilon_hb.template modify<LMPHostType>();
  k_a_hb.template modify<LMPHostType>();
  k_cut_hb_0.template modify<LMPHostType>();
  k_cut_hb_c.template modify<LMPHostType>();
  k_cut_hb_lo.template modify<LMPHostType>();
  k_cut_hb_hi.template modify<LMPHostType>();
  k_cut_hb_lc.template modify<LMPHostType>();
  k_cut_hb_hc.template modify<LMPHostType>();
  k_b_hb_lo.template modify<LMPHostType>();
  k_b_hb_hi.template modify<LMPHostType>();
  k_shift_hb.template modify<LMPHostType>();
  k_cutsq_hb_hc.template modify<LMPHostType>();

  k_a_hb1.template modify<LMPHostType>();
  k_theta_hb1_0.template modify<LMPHostType>();
  k_dtheta_hb1_ast.template modify<LMPHostType>();
  k_b_hb1.template modify<LMPHostType>();
  k_dtheta_hb1_c.template modify<LMPHostType>();

  k_a_hb2.template modify<LMPHostType>();
  k_theta_hb2_0.template modify<LMPHostType>();
  k_dtheta_hb2_ast.template modify<LMPHostType>();
  k_b_hb2.template modify<LMPHostType>();
  k_dtheta_hb2_c.template modify<LMPHostType>();

  k_a_hb3.template modify<LMPHostType>();
  k_theta_hb3_0.template modify<LMPHostType>();
  k_dtheta_hb3_ast.template modify<LMPHostType>();
  k_b_hb3.template modify<LMPHostType>();
  k_dtheta_hb3_c.template modify<LMPHostType>();

  k_a_hb4.template modify<LMPHostType>();
  k_theta_hb4_0.template modify<LMPHostType>();
  k_dtheta_hb4_ast.template modify<LMPHostType>();
  k_b_hb4.template modify<LMPHostType>();
  k_dtheta_hb4_c.template modify<LMPHostType>();

  k_a_hb7.template modify<LMPHostType>();
  k_theta_hb7_0.template modify<LMPHostType>();
  k_dtheta_hb7_ast.template modify<LMPHostType>();
  k_b_hb7.template modify<LMPHostType>();
  k_dtheta_hb7_c.template modify<LMPHostType>();

  k_a_hb8.template modify<LMPHostType>();
  k_theta_hb8_0.template modify<LMPHostType>();
  k_dtheta_hb8_ast.template modify<LMPHostType>();
  k_b_hb8.template modify<LMPHostType>();
  k_dtheta_hb8_c.template modify<LMPHostType>();

  // Sync to device
  k_epsilon_hb.template sync<DeviceType>();
  k_a_hb.template sync<DeviceType>();
  k_cut_hb_0.template sync<DeviceType>();
  k_cut_hb_c.template sync<DeviceType>();
  k_cut_hb_lo.template sync<DeviceType>();
  k_cut_hb_hi.template sync<DeviceType>();
  k_cut_hb_lc.template sync<DeviceType>();
  k_cut_hb_hc.template sync<DeviceType>();
  k_b_hb_lo.template sync<DeviceType>();
  k_b_hb_hi.template sync<DeviceType>();
  k_shift_hb.template sync<DeviceType>();
  k_cutsq_hb_hc.template sync<DeviceType>();

  k_a_hb1.template sync<DeviceType>();
  k_theta_hb1_0.template sync<DeviceType>();
  k_dtheta_hb1_ast.template sync<DeviceType>();
  k_b_hb1.template sync<DeviceType>();
  k_dtheta_hb1_c.template sync<DeviceType>();

  k_a_hb2.template sync<DeviceType>();
  k_theta_hb2_0.template sync<DeviceType>();
  k_dtheta_hb2_ast.template sync<DeviceType>();
  k_b_hb2.template sync<DeviceType>();
  k_dtheta_hb2_c.template sync<DeviceType>();

  k_a_hb3.template sync<DeviceType>();
  k_theta_hb3_0.template sync<DeviceType>();
  k_dtheta_hb3_ast.template sync<DeviceType>();
  k_b_hb3.template sync<DeviceType>();
  k_dtheta_hb3_c.template sync<DeviceType>();

  k_a_hb4.template sync<DeviceType>();
  k_theta_hb4_0.template sync<DeviceType>();
  k_dtheta_hb4_ast.template sync<DeviceType>();
  k_b_hb4.template sync<DeviceType>();
  k_dtheta_hb4_c.template sync<DeviceType>();

  k_a_hb7.template sync<DeviceType>();
  k_theta_hb7_0.template sync<DeviceType>();
  k_dtheta_hb7_ast.template sync<DeviceType>();
  k_b_hb7.template sync<DeviceType>();
  k_dtheta_hb7_c.template sync<DeviceType>();

  k_a_hb8.template sync<DeviceType>();
  k_theta_hb8_0.template sync<DeviceType>();
  k_dtheta_hb8_ast.template sync<DeviceType>();
  k_b_hb8.template sync<DeviceType>();
  k_dtheta_hb8_c.template sync<DeviceType>();

  // Register the COM screen cutoff for this pair: the h-bond cutoff acts at the
  // base site (COM +/- 0.4*nx on each atom), so a COM-distance screen needs a
  // 2*0.4 margin to never drop an interacting pair. The npair fix takes the max
  // over all consuming styles and type-pairs.
  if (fix_oxdna_npairKK) fix_oxdna_npairKK->request_screen_cutoff(cutone + 0.8);

  // "cutone" is "cut_hb_hc[i][j]", sets the master list distance cutoff
  return cutone;

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
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
int PairOxdnaHbondKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairOxdnaHbondKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdnaHbondKokkos<LMPHostType>;
#endif
}
