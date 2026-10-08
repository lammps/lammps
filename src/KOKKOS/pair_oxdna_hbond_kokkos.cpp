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
#include "pair_oxdna_hbond_kokkos_impl.h"
#include "pair_oxdna3_xstk_kokkos.h"

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
  screened_launch_count = 0;
  fuse_partner = nullptr;
  fuse_ncompute = 0;
  unique_basepair_enabled = 0;
  last_idc_nbuild = -1;
  last_idc_nall = -1;
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

  x = fix_oxdna_lrfKK->packed_x();
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

  if (execution_space == HostKK) {
    NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
    d_neighbors = k_list->d_neighbors;
    anum = list->inum;
    d_alist = k_list->d_ilist;
    d_numneigh = k_list->d_numneigh;
  }

  int need_dup = lmp->kokkos->need_dup<DeviceType>();
  if (need_dup) {
    fix_oxdna_lrfKK->prepare_dup_f_torque(f, torque);
    dup_f = fix_oxdna_lrfKK->dup_f;
    dup_torque = fix_oxdna_lrfKK->dup_torque;
  } else {
    ndup_f = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_torque = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, \
    Kokkos::Experimental::ScatterNonDuplicated>(torque);
  }
  if (eflag_atom) {
    if (need_dup)
      dup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    else
      ndup_eatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
  }
  if (vflag_atom) {
    if (need_dup)
      dup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterDuplicated>(d_vatom);
    else
      ndup_vatom = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum,
        Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  copymode = 1;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->packed_nx();
  d_ny_xtrct = fix_oxdna_lrfKK->packed_ny();
  d_nz_xtrct = fix_oxdna_lrfKK->packed_nz();

  // If we're on a GPU, look up fix_oxdna_npairKK screened pair count and pair_a/b views.
  if (execution_space != HostKK) {
    screened_pair_count = fix_oxdna_npairKK->screened_pair_count;
    d_pairs_screened = fix_oxdna_npairKK->k_pairs_screened.template view<DeviceType>();
    d_screened_offsets = fix_oxdna_npairKK->get_screened_offsets();
#if OXDNA_KK_SCREENED_PER_ATOM
    screened_launch_count = fix_oxdna_npairKK->get_anum();
#else
    screened_launch_count = screened_pair_count;
#endif
#if OXDNA_KK_TWO_PHASE
    if (static_cast<int>(d_radial_pairs.extent(0)) < screened_pair_count)
      d_radial_pairs = typename AT::t_int_1d(Kokkos::view_alloc(Kokkos::WithoutInitializing, "pair:radial_pairs"),
                                             screened_pair_count + screened_pair_count/10);
    if (d_radial_count.data() == nullptr) d_radial_count = typename AT::t_int_scalar("pair:radial_count");
    Kokkos::deep_copy(d_radial_count, 0);
#endif
  }

  // the complementary nucleotide IDs are a custom per-atom vector that is
  // reallocated when the per-atom arrays grow and whose ghost values are only
  // updated by border communication, so re-fetch the pointer every step and
  // only copy the values when the neighbor list was rebuilt

  idc = unique_basepair_enabled ? atom->ivector[idc_index] : nullptr;
  if (idc) {
    const int nall = atom->nlocal + atom->nghost;
    if ((neighbor->nbuild != last_idc_nbuild) || (nall != last_idc_nall)) {
      if (k_idc.extent(0) < static_cast<size_t>(nall))
        k_idc = DAT::tdual_int_1d("pair:idc", nall);
      atomKK->sync(Host, IVECTOR_MASK);
      auto h_idc = k_idc.view_host();
      for (int i = 0; i < nall; i++) h_idc(i) = idc[i];
      k_idc.modify_host();
      k_idc.template sync<DeviceType>();
      last_idc_nbuild = neighbor->nbuild;
      last_idc_nall = nall;
    }
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
  auto run_compute_gpu = [&](auto gpu_tag, auto radial_tag, auto evflag_tag) {
    constexpr int EVFLAG = decltype(evflag_tag)::value;
#if OXDNA_KK_TWO_PHASE
    Kokkos::parallel_for(OxdnaPairRangePolicy<DeviceType, decltype(radial_tag)>(0,screened_pair_count),*this);
#else
    (void) radial_tag;
#endif
    if constexpr (EVFLAG) {
      Kokkos::parallel_reduce(OxdnaPairRangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_launch_count),*this,ev);
    } else {
      Kokkos::parallel_for(OxdnaPairRangePolicy<DeviceType, decltype(gpu_tag)>(0,screened_launch_count),*this);
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
        run_compute_gpu(TagPairOxdnaHbondComputeGPUPair<OXDNA,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{},
                        TagPairOxdnaHbondComputeGPURadial<OXDNA,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      }
    } else if (oxdnaflag == OXDNA3) {
      if (use_host_launch) {
        run_compute_host(TagPairOxdnaHbondCompute<OXDNA3,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
      } else {
        run_compute_gpu(TagPairOxdnaHbondComputeGPUPair<OXDNA3,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{},
                        TagPairOxdnaHbondComputeGPURadial<OXDNA3,NEIGHFLAG,NEWTON_PAIR,EVFLAG>{}, evflag_tag);
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

  // with the fused hbond + oxdna3/xstk kernel, the style that is computed
  // second in this force evaluation launches it for both styles
  const bool fused = (fuse_partner != nullptr) && !eflag_atom && !vflag_atom && !need_dup;
  if (fused) {
    fuse_ncompute++;
    if (fuse_partner->fuse_ncompute == fuse_ncompute)
      fuse_partner->compute_fused(this, fuse_partner);
    else if (fuse_partner->fuse_ncompute != fuse_ncompute - 1)
      error->one(FLERR, "Fused kernel of pair oxdna3/hbond/kk and oxdna3/xstk/kk out of step");
  } else {
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
  }

  if (need_dup) {
    Kokkos::Experimental::contribute(f, dup_f);
    Kokkos::Experimental::contribute(torque, dup_torque);
  }

  if (eflag_global) eng_vdwl += static_cast<double>(ev.evdwl);
  if (vflag_global) {
    virial[0] += static_cast<double>(ev.v[0]);
    virial[1] += static_cast<double>(ev.v[1]);
    virial[2] += static_cast<double>(ev.v[2]);
    virial[3] += static_cast<double>(ev.v[3]);
    virial[4] += static_cast<double>(ev.v[4]);
    virial[5] += static_cast<double>(ev.v[5]);
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

  if (eflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_eatom, dup_eatom);
    k_eatom.template modify<DeviceType>();
    k_eatom.sync_host();
  }

  if (vflag_atom) {
    if (need_dup)
      Kokkos::Experimental::contribute(d_vatom, dup_vatom);
    k_vatom.template modify<DeviceType>();
    k_vatom.sync_host();
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
    constexpr KK_FLOAT dx_cbs_oxdna1=static_cast<KK_FLOAT>(+0.4);
    ra_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
  } else if (OXDNAFLAG == OXDNA3) {
    constexpr KK_FLOAT dx_cbs_pur_oxdna3 = static_cast<KK_FLOAT>(+0.43);
    constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = static_cast<KK_FLOAT>(+0.37);
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
    const KK_FLOAT factor_lj = static_cast<KK_FLOAT>(special_lj[sbmask(b)]);
    b &= NEIGHMASK;
    const int btype = type(b);

    // no hydrogen bonding between these base types (e.g. non-complementary bases):
    // f1 and thus the energy would be zero, so skip the site geometry altogether
    if (d_params_hb(atype,btype).epsilon_hb == static_cast<KK_FLOAT>(0.0)) continue;

    if (unique_basepair_enabled) {
      const int idca = d_idc(a);
      const int idcb = d_idc(b);
      if (idca != tag(b) && idcb != tag(a) && idca > 0 && idcb > 0) continue;
    }

    // vector COM-hbond site b
    if (OXDNAFLAG == OXDNA) {
      // Used for oxdna1 and oxdna2, which have the same hbond site offset.
      constexpr KK_FLOAT dx_cbs_oxdna1=static_cast<KK_FLOAT>(+0.4);
      rb_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(b,0);
      rb_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(b,1);
      rb_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(b,2);
    } else if (OXDNAFLAG == OXDNA3) {
      constexpr KK_FLOAT dx_cbs_pur_oxdna3 = static_cast<KK_FLOAT>(+0.43);
      constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = static_cast<KK_FLOAT>(+0.37);
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
    r_hb = Kokkos::sqrt(rsq_hb);
    rinv_hb = static_cast<KK_FLOAT>(1.0) / r_hb;

    delr_hb_norm[0] = delr_hb[0] * rinv_hb;
    delr_hb_norm[1] = delr_hb[1] * rinv_hb;
    delr_hb_norm[2] = delr_hb[2] * rinv_hb;

    // beginning of modulation factors

    // f1 = f1 modulation factor
    f1 = F1_KK(r_hb, d_params_hb(atype,btype).epsilon_hb, d_params_hb(atype,btype).a_hb, d_params_hb(atype,btype).cut_hb_0,
            d_params_hb(atype,btype).cut_hb_lc, d_params_hb(atype,btype).cut_hb_hc, d_params_hb(atype,btype).cut_hb_lo,
            d_params_hb(atype,btype).cut_hb_hi, d_params_hb(atype,btype).b_hb_lo,
            d_params_hb(atype,btype).b_hb_hi, d_params_hb(atype,btype).shift_hb);

    // start early rejection criterium
    if (f1) {
      // theta1 calculation
      cost1 = - (d_nx_xtrct(a,0)*d_nx_xtrct(b,0) + d_nx_xtrct(a,1)*d_nx_xtrct(b,1) + d_nx_xtrct(a,2)*d_nx_xtrct(b,2));
      if (cost1 > static_cast<KK_FLOAT>(1.0)) cost1 = static_cast<KK_FLOAT>(1.0);
      if (cost1 < static_cast<KK_FLOAT>(-1.0)) cost1 = static_cast<KK_FLOAT>(-1.0);
      theta1 = Kokkos::acos(cost1);
      // f4t1 = f4 modulation factor
      f4t1 = F4_KK(theta1, d_params_hb(atype,btype).a_hb1, d_params_hb(atype, btype).theta_hb1_0, d_params_hb(atype, btype).dtheta_hb1_ast,
              d_params_hb(atype, btype).b_hb1, d_params_hb(atype, btype).dtheta_hb1_c);
    // end of f1

    // f4t1 early rejection criterium
    if (f4t1) {
      // theta2 calculation
      cost2 = - (d_nx_xtrct(a,0)*delr_hb_norm[0] + d_nx_xtrct(a,1)*delr_hb_norm[1] + d_nx_xtrct(a,2)*delr_hb_norm[2]);
      if (cost2 > static_cast<KK_FLOAT>(1.0)) cost2 = static_cast<KK_FLOAT>(1.0);
      if (cost2 < static_cast<KK_FLOAT>(-1.0)) cost2 = static_cast<KK_FLOAT>(-1.0);
      theta2 = Kokkos::acos(cost2);
      // f4t2 = f4 modulation factor
      f4t2 = F4_KK(theta2, d_params_hb(atype,btype).a_hb2, d_params_hb(atype, btype).theta_hb2_0, d_params_hb(atype, btype).dtheta_hb2_ast,
              d_params_hb(atype, btype).b_hb2, d_params_hb(atype, btype).dtheta_hb2_c);
    // end of f4t1

    // f4t2 early rejection criterium
    if (f4t2) {
      // theta3 calculation
      cost3 = d_nx_xtrct(b,0)*delr_hb_norm[0] + d_nx_xtrct(b,1)*delr_hb_norm[1] + d_nx_xtrct(b,2)*delr_hb_norm[2];
      if (cost3 > static_cast<KK_FLOAT>(1.0)) cost3 = static_cast<KK_FLOAT>(1.0);
      if (cost3 < static_cast<KK_FLOAT>(-1.0)) cost3 = static_cast<KK_FLOAT>(-1.0);
      theta3 = Kokkos::acos(cost3);
      // f4t3 = f4 modulation factor
      f4t3 = F4_KK(theta3, d_params_hb(atype,btype).a_hb3, d_params_hb(atype, btype).theta_hb3_0, d_params_hb(atype, btype).dtheta_hb3_ast,
              d_params_hb(atype, btype).b_hb3, d_params_hb(atype, btype).dtheta_hb3_c);
    // end of f4t2

    // f4t3 early rejection criterium
    if (f4t3) {
      // theta4 calculation
      cost4 = d_nz_xtrct(a,0)*d_nz_xtrct(b,0) + d_nz_xtrct(a,1)*d_nz_xtrct(b,1) + d_nz_xtrct(a,2)*d_nz_xtrct(b,2);
      if (cost4 > static_cast<KK_FLOAT>(1.0)) cost4 = static_cast<KK_FLOAT>(1.0);
      if (cost4 < static_cast<KK_FLOAT>(-1.0)) cost4 = static_cast<KK_FLOAT>(-1.0);
      theta4 = Kokkos::acos(cost4);
      // f4t4 = f4 modulation factor
      f4t4 = F4_KK(theta4, d_params_hb(atype,btype).a_hb4, d_params_hb(atype, btype).theta_hb4_0, d_params_hb(atype, btype).dtheta_hb4_ast,
              d_params_hb(atype, btype).b_hb4, d_params_hb(atype, btype).dtheta_hb4_c);
    // end of f4t3

    // f4t4 early rejection criterium
    if (f4t4) {
      cost7 = - (d_nz_xtrct(a,0)*delr_hb_norm[0] + d_nz_xtrct(a,1)*delr_hb_norm[1] + d_nz_xtrct(a,2)*delr_hb_norm[2]);
      if (cost7 > static_cast<KK_FLOAT>(1.0)) cost7 = static_cast<KK_FLOAT>(1.0);
      if (cost7 < static_cast<KK_FLOAT>(-1.0)) cost7 = static_cast<KK_FLOAT>(-1.0);
      theta7 = Kokkos::acos(cost7);
      // f4t7 = f4 modulation factor
      f4t7 = F4_KK(theta7, d_params_hb(atype,btype).a_hb7, d_params_hb(atype, btype).theta_hb7_0, d_params_hb(atype, btype).dtheta_hb7_ast,
              d_params_hb(atype, btype).b_hb7, d_params_hb(atype, btype).dtheta_hb7_c);
    // end of f4t4

    // f4t7 early rejection criterium
    if (f4t7) {
      cost8 = d_nz_xtrct(b,0)*delr_hb_norm[0] + d_nz_xtrct(b,1)*delr_hb_norm[1] + d_nz_xtrct(b,2)*delr_hb_norm[2];
      if (cost8 > static_cast<KK_FLOAT>(1.0)) cost8 = static_cast<KK_FLOAT>(1.0);
      if (cost8 < static_cast<KK_FLOAT>(-1.0)) cost8 = static_cast<KK_FLOAT>(-1.0);
      theta8 = Kokkos::acos(cost8);
      // f4t8 = f4 modulation factor
      f4t8 = F4_KK(theta8, d_params_hb(atype,btype).a_hb8, d_params_hb(atype, btype).theta_hb8_0, d_params_hb(atype, btype).dtheta_hb8_ast,
              d_params_hb(atype, btype).b_hb8, d_params_hb(atype, btype).dtheta_hb8_c);

      evdwl = static_cast<KK_ACC_FLOAT>(f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj);
    // end of f4t7

    // evdwl early rejection criterium
    if (evdwl) {
      // df1 = DF1 modulation factor
      df1 = DF1_KK(r_hb, d_params_hb(atype,btype).epsilon_hb, d_params_hb(atype,btype).a_hb, d_params_hb(atype,btype).cut_hb_0,
            d_params_hb(atype,btype).cut_hb_lc, d_params_hb(atype,btype).cut_hb_hc, d_params_hb(atype,btype).cut_hb_lo,
            d_params_hb(atype,btype).cut_hb_hi, d_params_hb(atype,btype).b_hb_lo,
            d_params_hb(atype,btype).b_hb_hi);
      // df4t1 = DF4 modulation factor
      df4t1 = DF4_KK(theta1, d_params_hb(atype,btype).a_hb1, d_params_hb(atype, btype).theta_hb1_0, d_params_hb(atype, btype).dtheta_hb1_ast,
              d_params_hb(atype, btype).b_hb1, d_params_hb(atype, btype).dtheta_hb1_c)/Kokkos::sin(theta1);
      // df4t2 = DF4 modulation factor
      df4t2 = DF4_KK(theta2, d_params_hb(atype,btype).a_hb2, d_params_hb(atype, btype).theta_hb2_0, d_params_hb(atype, btype).dtheta_hb2_ast,
              d_params_hb(atype, btype).b_hb2, d_params_hb(atype, btype).dtheta_hb2_c)/Kokkos::sin(theta2);
      // df4t3 = DF4 modulation factor
      df4t3 = DF4_KK(theta3, d_params_hb(atype,btype).a_hb3, d_params_hb(atype, btype).theta_hb3_0, d_params_hb(atype, btype).dtheta_hb3_ast,
              d_params_hb(atype, btype).b_hb3, d_params_hb(atype, btype).dtheta_hb3_c)/Kokkos::sin(theta3);
      // df4t4 = DF4 modulation factor
      df4t4 = DF4_KK(theta4, d_params_hb(atype,btype).a_hb4, d_params_hb(atype, btype).theta_hb4_0, d_params_hb(atype, btype).dtheta_hb4_ast,
              d_params_hb(atype, btype).b_hb4, d_params_hb(atype, btype).dtheta_hb4_c)/Kokkos::sin(theta4);
      // df4t7 = DF4 modulation factor
      df4t7 = DF4_KK(theta7, d_params_hb(atype,btype).a_hb7, d_params_hb(atype, btype).theta_hb7_0, d_params_hb(atype, btype).dtheta_hb7_ast,
              d_params_hb(atype, btype).b_hb7, d_params_hb(atype, btype).dtheta_hb7_c)/Kokkos::sin(theta7);
      // df4t8 = DF4 modulation factor
      df4t8 = DF4_KK(theta8, d_params_hb(atype,btype).a_hb8, d_params_hb(atype, btype).theta_hb8_0, d_params_hb(atype, btype).dtheta_hb8_ast,
              d_params_hb(atype, btype).b_hb8, d_params_hb(atype, btype).dtheta_hb8_c)/Kokkos::sin(theta8);

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
      finc  = static_cast<KK_ACC_FLOAT>(-df1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj);

      delf[0] += static_cast<KK_ACC_FLOAT>(delr_hb[0]) * finc;
      delf[1] += static_cast<KK_ACC_FLOAT>(delr_hb[1]) * finc;
      delf[2] += static_cast<KK_ACC_FLOAT>(delr_hb[2]) * finc;

      // theta2 force
      if (theta2 != static_cast<KK_FLOAT>(0.0)) {

        finc  = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj);

        delf[0] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[0]*cost2 + d_nx_xtrct(a,0)) * finc;
        delf[1] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[1]*cost2 + d_nx_xtrct(a,1)) * finc;
        delf[2] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[2]*cost2 + d_nx_xtrct(a,2)) * finc;
      }

      // theta3 force
      if (theta3 != static_cast<KK_FLOAT>(0.0)) {

        finc  = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * rinv_hb * factor_lj);

        delf[0] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[0]*cost3 - d_nx_xtrct(b,0)) * finc;
        delf[1] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[1]*cost3 - d_nx_xtrct(b,1)) * finc;
        delf[2] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[2]*cost3 - d_nx_xtrct(b,2)) * finc;
      }

      // theta7 force
      if (theta7 != static_cast<KK_FLOAT>(0.0)) {

        finc  = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * rinv_hb * factor_lj);

        delf[0] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[0]*cost7 + d_nz_xtrct(a,0)) * finc;
        delf[1] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[1]*cost7 + d_nz_xtrct(a,1)) * finc;
        delf[2] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[2]*cost7 + d_nz_xtrct(a,2)) * finc;

      }

      // theta8 force
      if (theta8 != static_cast<KK_FLOAT>(0.0)) {

        finc  = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * rinv_hb * factor_lj);

        delf[0] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[0]*cost8 - d_nz_xtrct(b,0)) * finc;
        delf[1] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[1]*cost8 - d_nz_xtrct(b,1)) * finc;
        delf[2] += static_cast<KK_ACC_FLOAT>(delr_hb_norm[2]*cost8 - d_nz_xtrct(b,2)) * finc;

      }

      // increment forces and torques

      a_f(a,0) += delf[0];
      a_f(a,1) += delf[1];
      a_f(a,2) += delf[2];
      delta[0] = static_cast<KK_ACC_FLOAT>(ra_chb[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_chb[2])*delf[1];
      delta[1] = static_cast<KK_ACC_FLOAT>(ra_chb[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_chb[0])*delf[2];
      delta[2] = static_cast<KK_ACC_FLOAT>(ra_chb[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_chb[1])*delf[0];
      a_torque(a,0) += delta[0];
      a_torque(a,1) += delta[1];
      a_torque(a,2) += delta[2];

      if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
        a_f(b,0) -= delf[0];
        a_f(b,1) -= delf[1];
        a_f(b,2) -= delf[2];
        deltb[0] = static_cast<KK_ACC_FLOAT>(rb_chb[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_chb[2])*delf[1];
        deltb[1] = static_cast<KK_ACC_FLOAT>(rb_chb[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_chb[0])*delf[2];
        deltb[2] = static_cast<KK_ACC_FLOAT>(rb_chb[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_chb[1])*delf[0];
        a_torque(b,0) -= deltb[0];
        a_torque(b,1) -= deltb[1];
        a_torque(b,2) -= deltb[2];
      }

      // increment energy and virial
      // NOTE: The virial is calculated on the 'molecular' basis.
      // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

      if (EVFLAG) {
        ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?static_cast<KK_ACC_FLOAT>(1.0):static_cast<KK_ACC_FLOAT>(0.5))*evdwl;

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
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

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * df4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj);

        t1dir[0] = d_nx_xtrct(a,1) * d_nx_xtrct(b,2) - d_nx_xtrct(a,2) * d_nx_xtrct(b,1);
        t1dir[1] = d_nx_xtrct(a,2) * d_nx_xtrct(b,0) - d_nx_xtrct(a,0) * d_nx_xtrct(b,2);
        t1dir[2] = d_nx_xtrct(a,0) * d_nx_xtrct(b,1) - d_nx_xtrct(a,1) * d_nx_xtrct(b,0);
        delta[0] += static_cast<KK_ACC_FLOAT>(t1dir[0]) * tpair;
        delta[1] += static_cast<KK_ACC_FLOAT>(t1dir[1]) * tpair;
        delta[2] += static_cast<KK_ACC_FLOAT>(t1dir[2]) * tpair;
        deltb[0] += static_cast<KK_ACC_FLOAT>(t1dir[0]) * tpair;
        deltb[1] += static_cast<KK_ACC_FLOAT>(t1dir[1]) * tpair;
        deltb[2] += static_cast<KK_ACC_FLOAT>(t1dir[2]) * tpair;
      }
      //theta2 torque
      if (theta2 != static_cast<KK_FLOAT>(0.0)) {

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * df4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj);

        t2dir[0] = d_nx_xtrct(a,1) * delr_hb_norm[2] - d_nx_xtrct(a,2) * delr_hb_norm[1];
        t2dir[1] = d_nx_xtrct(a,2) * delr_hb_norm[0] - d_nx_xtrct(a,0) * delr_hb_norm[2];
        t2dir[2] = d_nx_xtrct(a,0) * delr_hb_norm[1] - d_nx_xtrct(a,1) * delr_hb_norm[0];
        delta[0] += static_cast<KK_ACC_FLOAT>(t2dir[0]) * tpair;
        delta[1] += static_cast<KK_ACC_FLOAT>(t2dir[1]) * tpair;
        delta[2] += static_cast<KK_ACC_FLOAT>(t2dir[2]) * tpair;
      }
      //theta3 torque
      if (theta3 != static_cast<KK_FLOAT>(0.0)) {

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * df4t3 * f4t4 * f4t7 * f4t8 * factor_lj);

        t3dir[0] = d_nx_xtrct(b,1) * delr_hb_norm[2] - d_nx_xtrct(b,2) * delr_hb_norm[1];
        t3dir[1] = d_nx_xtrct(b,2) * delr_hb_norm[0] - d_nx_xtrct(b,0) * delr_hb_norm[2];
        t3dir[2] = d_nx_xtrct(b,0) * delr_hb_norm[1] - d_nx_xtrct(b,1) * delr_hb_norm[0];
        deltb[0] += static_cast<KK_ACC_FLOAT>(t3dir[0]) * tpair;
        deltb[1] += static_cast<KK_ACC_FLOAT>(t3dir[1]) * tpair;
        deltb[2] += static_cast<KK_ACC_FLOAT>(t3dir[2]) * tpair;
      }
      //theta4 torque
      if (theta4 != static_cast<KK_FLOAT>(0.0)) {

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * f4t3 * df4t4 * f4t7 * f4t8 * factor_lj);

        t4dir[0] = d_nz_xtrct(b,1) * d_nz_xtrct(a,2) - d_nz_xtrct(b,2) * d_nz_xtrct(a,1);
        t4dir[1] = d_nz_xtrct(b,2) * d_nz_xtrct(a,0) - d_nz_xtrct(b,0) * d_nz_xtrct(a,2);
        t4dir[2] = d_nz_xtrct(b,0) * d_nz_xtrct(a,1) - d_nz_xtrct(b,1) * d_nz_xtrct(a,0);
        delta[0] += static_cast<KK_ACC_FLOAT>(t4dir[0]) * tpair;
        delta[1] += static_cast<KK_ACC_FLOAT>(t4dir[1]) * tpair;
        delta[2] += static_cast<KK_ACC_FLOAT>(t4dir[2]) * tpair;
        deltb[0] += static_cast<KK_ACC_FLOAT>(t4dir[0]) * tpair;
        deltb[1] += static_cast<KK_ACC_FLOAT>(t4dir[1]) * tpair;
        deltb[2] += static_cast<KK_ACC_FLOAT>(t4dir[2]) * tpair;
      }
      //theta7 torque
      if (theta7 != static_cast<KK_FLOAT>(0.0)) {

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * f4t3 * f4t4 * df4t7 * f4t8 * factor_lj);

        t7dir[0] = d_nz_xtrct(a,1) * delr_hb_norm[2] - d_nz_xtrct(a,2) * delr_hb_norm[1];
        t7dir[1] = d_nz_xtrct(a,2) * delr_hb_norm[0] - d_nz_xtrct(a,0) * delr_hb_norm[2];
        t7dir[2] = d_nz_xtrct(a,0) * delr_hb_norm[1] - d_nz_xtrct(a,1) * delr_hb_norm[0];
        delta[0] += static_cast<KK_ACC_FLOAT>(t7dir[0]) * tpair;
        delta[1] += static_cast<KK_ACC_FLOAT>(t7dir[1]) * tpair;
        delta[2] += static_cast<KK_ACC_FLOAT>(t7dir[2]) * tpair;
      }
      //theta8 torque
      if (theta8 != static_cast<KK_FLOAT>(0.0)) {

        tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * df4t8 * factor_lj);

        t8dir[0] = d_nz_xtrct(b,1) * delr_hb_norm[2] - d_nz_xtrct(b,2) * delr_hb_norm[1];
        t8dir[1] = d_nz_xtrct(b,2) * delr_hb_norm[0] - d_nz_xtrct(b,0) * delr_hb_norm[2];
        t8dir[2] = d_nz_xtrct(b,0) * delr_hb_norm[1] - d_nz_xtrct(b,1) * delr_hb_norm[0];
        deltb[0] += static_cast<KK_ACC_FLOAT>(t8dir[0]) * tpair;
        deltb[1] += static_cast<KK_ACC_FLOAT>(t8dir[1]) * tpair;
        deltb[2] += static_cast<KK_ACC_FLOAT>(t8dir[2]) * tpair;
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
   first phase of the two-phase evaluation: append the screened pairs that
   pass the radial test to d_radial_pairs
------------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondComputeGPURadial<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>,
  const int &ipair) const
{
  KK_ACC_FLOAT fa[3], ta[3];
  EV_FLOAT ev;
  if (this->template screened_pair_body<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG,1>(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), ipair, fa, ta, ev))
    d_radial_pairs(Kokkos::atomic_fetch_add(&d_radial_count(), 1)) = ipair;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaHbondKokkos<DeviceType>::operator()(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ipair, EV_FLOAT &ev) const
{
  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<Kokkos::Experimental::ScatterAtomic>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<Kokkos::Experimental::ScatterAtomic>();

  KK_ACC_FLOAT fa[3] = {0.0, 0.0, 0.0}, ta[3] = {0.0, 0.0, 0.0};
#if OXDNA_KK_SCREENED_PER_ATOM
  // one thread per atom: loop over its screened pairs, index ipair is the atom
  const int ibeg = d_screened_offsets(ipair);
  const int iend = d_screened_offsets(ipair+1);
  bool any = false;
  for (int jpair = ibeg; jpair < iend; jpair++)
    if (screened_pair_body(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), jpair, fa, ta, ev)) any = true;
  if (any) {
    const int a = static_cast<int>(d_pairs_screened(ibeg) >> 32);
    a_f(a,0) += fa[0];
    a_f(a,1) += fa[1];
    a_f(a,2) += fa[2];
    a_torque(a,0) += ta[0];
    a_torque(a,1) += ta[1];
    a_torque(a,2) += ta[2];
  }
#else
#if OXDNA_KK_TWO_PHASE
  // one thread per screened pair that passed the radial test
  if (ipair >= d_radial_count()) return;
  const int jpair = d_radial_pairs(ipair);
#else
  // one thread per screened pair
  const int jpair = ipair;
#endif
  if (screened_pair_body(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(), jpair, fa, ta, ev)) {
    const int a = static_cast<int>(d_pairs_screened(jpair) >> 32);
    a_f(a,0) += fa[0];
    a_f(a,1) += fa[1];
    a_f(a,2) += fa[2];
    a_torque(a,0) += ta[0];
    a_torque(a,1) += ta[1];
    a_torque(a,2) += ta[2];
  }
#endif
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

  k_params_hb = decltype(k_params_hb)("PairOxdnaHbondKokkos:params_hb", n+1, n+1);







  d_params_hb = k_params_hb.template view<DeviceType>();







}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR, "The oxDNA and oxRNA pair styles do not take any arguments");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaHbondKokkos<DeviceType>::init_style()
{
  // the internal helper fixes are always created for the default KOKKOS variant,
  // so /kk/host styles cannot work with them when LAMMPS is compiled for a GPU

  if (std::is_same_v<DeviceType, LMPHostType> && !std::is_same_v<DeviceType, LMPDeviceType>)
    error->all(FLERR, "The /kk/host variants of the CG-DNA styles are not supported "
               "when LAMMPS is compiled for a GPU");

  // on GPUs the screened-pair kernel runs over the pair list of fix
  // OXDNA/NPAIR/kk, so only the host kernel needs a neighbor list

  neighflag = lmp->kokkos->neighflag;
  if (execution_space == HostKK) {
    neighbor->add_request(this);
    auto request = neighbor->find_request(this);
    request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                             !std::is_same_v<DeviceType,LMPDeviceType>);
    request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
    if (neighflag == FULL) request->enable_full();
  }

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
  idc_index = -1;
  last_idc_nbuild = -1;
  last_idc_nall = -1;
  const int ifix = modify->find_fix("Basepairs");
  if (ifix >= 0) {
    int idx, flag, cols;
    idx = atom->find_custom("idc", flag, cols);
    if (idx >= 0 && flag == 0) {
      idc_index = idx;
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
  k_params_hb.view_host()(i,j).epsilon_hb = static_cast<KK_FLOAT>(epsilon_hb[i][j]); k_params_hb.view_host()(j,i).epsilon_hb = static_cast<KK_FLOAT>(epsilon_hb[j][i]);
  k_params_hb.view_host()(i,j).a_hb = static_cast<KK_FLOAT>(a_hb[i][j]); k_params_hb.view_host()(j,i).a_hb = static_cast<KK_FLOAT>(a_hb[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_0 = static_cast<KK_FLOAT>(cut_hb_0[i][j]); k_params_hb.view_host()(j,i).cut_hb_0 = static_cast<KK_FLOAT>(cut_hb_0[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_c = static_cast<KK_FLOAT>(cut_hb_c[i][j]); k_params_hb.view_host()(j,i).cut_hb_c = static_cast<KK_FLOAT>(cut_hb_c[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_lo = static_cast<KK_FLOAT>(cut_hb_lo[i][j]); k_params_hb.view_host()(j,i).cut_hb_lo = static_cast<KK_FLOAT>(cut_hb_lo[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_hi = static_cast<KK_FLOAT>(cut_hb_hi[i][j]); k_params_hb.view_host()(j,i).cut_hb_hi = static_cast<KK_FLOAT>(cut_hb_hi[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_lc = static_cast<KK_FLOAT>(cut_hb_lc[i][j]); k_params_hb.view_host()(j,i).cut_hb_lc = static_cast<KK_FLOAT>(cut_hb_lc[j][i]);
  k_params_hb.view_host()(i,j).cut_hb_hc = static_cast<KK_FLOAT>(cut_hb_hc[i][j]); k_params_hb.view_host()(j,i).cut_hb_hc = static_cast<KK_FLOAT>(cut_hb_hc[j][i]);
  k_params_hb.view_host()(i,j).b_hb_lo = static_cast<KK_FLOAT>(b_hb_lo[i][j]); k_params_hb.view_host()(j,i).b_hb_lo = static_cast<KK_FLOAT>(b_hb_lo[j][i]);
  k_params_hb.view_host()(i,j).b_hb_hi = static_cast<KK_FLOAT>(b_hb_hi[i][j]); k_params_hb.view_host()(j,i).b_hb_hi = static_cast<KK_FLOAT>(b_hb_hi[j][i]);
  k_params_hb.view_host()(i,j).shift_hb = static_cast<KK_FLOAT>(shift_hb[i][j]); k_params_hb.view_host()(j,i).shift_hb = static_cast<KK_FLOAT>(shift_hb[j][i]);
  k_params_hb.view_host()(i,j).cutsq_hb_hc = static_cast<KK_FLOAT>(cutsq_hb_hc[i][j]); k_params_hb.view_host()(j,i).cutsq_hb_hc = static_cast<KK_FLOAT>(cutsq_hb_hc[j][i]);

  k_params_hb.view_host()(i,j).a_hb1 = static_cast<KK_FLOAT>(a_hb1[i][j]); k_params_hb.view_host()(j,i).a_hb1 = static_cast<KK_FLOAT>(a_hb1[j][i]);
  k_params_hb.view_host()(i,j).theta_hb1_0 = static_cast<KK_FLOAT>(theta_hb1_0[i][j]); k_params_hb.view_host()(j,i).theta_hb1_0 = static_cast<KK_FLOAT>(theta_hb1_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb1_ast = static_cast<KK_FLOAT>(dtheta_hb1_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb1_ast = static_cast<KK_FLOAT>(dtheta_hb1_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb1 = static_cast<KK_FLOAT>(b_hb1[i][j]); k_params_hb.view_host()(j,i).b_hb1 = static_cast<KK_FLOAT>(b_hb1[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb1_c = static_cast<KK_FLOAT>(dtheta_hb1_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb1_c = static_cast<KK_FLOAT>(dtheta_hb1_c[j][i]);

  k_params_hb.view_host()(i,j).a_hb2 = static_cast<KK_FLOAT>(a_hb2[i][j]); k_params_hb.view_host()(j,i).a_hb2 = static_cast<KK_FLOAT>(a_hb2[j][i]);
  k_params_hb.view_host()(i,j).theta_hb2_0 = static_cast<KK_FLOAT>(theta_hb2_0[i][j]); k_params_hb.view_host()(j,i).theta_hb2_0 = static_cast<KK_FLOAT>(theta_hb2_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb2_ast = static_cast<KK_FLOAT>(dtheta_hb2_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb2_ast = static_cast<KK_FLOAT>(dtheta_hb2_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb2 = static_cast<KK_FLOAT>(b_hb2[i][j]); k_params_hb.view_host()(j,i).b_hb2 = static_cast<KK_FLOAT>(b_hb2[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb2_c = static_cast<KK_FLOAT>(dtheta_hb2_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb2_c = static_cast<KK_FLOAT>(dtheta_hb2_c[j][i]);

  k_params_hb.view_host()(i,j).a_hb3 = static_cast<KK_FLOAT>(a_hb3[i][j]); k_params_hb.view_host()(j,i).a_hb3 = static_cast<KK_FLOAT>(a_hb3[j][i]);
  k_params_hb.view_host()(i,j).theta_hb3_0 = static_cast<KK_FLOAT>(theta_hb3_0[i][j]); k_params_hb.view_host()(j,i).theta_hb3_0 = static_cast<KK_FLOAT>(theta_hb3_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb3_ast = static_cast<KK_FLOAT>(dtheta_hb3_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb3_ast = static_cast<KK_FLOAT>(dtheta_hb3_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb3 = static_cast<KK_FLOAT>(b_hb3[i][j]); k_params_hb.view_host()(j,i).b_hb3 = static_cast<KK_FLOAT>(b_hb3[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb3_c = static_cast<KK_FLOAT>(dtheta_hb3_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb3_c = static_cast<KK_FLOAT>(dtheta_hb3_c[j][i]);

  k_params_hb.view_host()(i,j).a_hb4 = static_cast<KK_FLOAT>(a_hb4[i][j]); k_params_hb.view_host()(j,i).a_hb4 = static_cast<KK_FLOAT>(a_hb4[j][i]);
  k_params_hb.view_host()(i,j).theta_hb4_0 = static_cast<KK_FLOAT>(theta_hb4_0[i][j]); k_params_hb.view_host()(j,i).theta_hb4_0 = static_cast<KK_FLOAT>(theta_hb4_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb4_ast = static_cast<KK_FLOAT>(dtheta_hb4_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb4_ast = static_cast<KK_FLOAT>(dtheta_hb4_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb4 = static_cast<KK_FLOAT>(b_hb4[i][j]); k_params_hb.view_host()(j,i).b_hb4 = static_cast<KK_FLOAT>(b_hb4[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb4_c = static_cast<KK_FLOAT>(dtheta_hb4_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb4_c = static_cast<KK_FLOAT>(dtheta_hb4_c[j][i]);

  k_params_hb.view_host()(i,j).a_hb7 = static_cast<KK_FLOAT>(a_hb7[i][j]); k_params_hb.view_host()(j,i).a_hb7 = static_cast<KK_FLOAT>(a_hb7[j][i]);
  k_params_hb.view_host()(i,j).theta_hb7_0 = static_cast<KK_FLOAT>(theta_hb7_0[i][j]); k_params_hb.view_host()(j,i).theta_hb7_0 = static_cast<KK_FLOAT>(theta_hb7_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb7_ast = static_cast<KK_FLOAT>(dtheta_hb7_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb7_ast = static_cast<KK_FLOAT>(dtheta_hb7_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb7 = static_cast<KK_FLOAT>(b_hb7[i][j]); k_params_hb.view_host()(j,i).b_hb7 = static_cast<KK_FLOAT>(b_hb7[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb7_c = static_cast<KK_FLOAT>(dtheta_hb7_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb7_c = static_cast<KK_FLOAT>(dtheta_hb7_c[j][i]);

  k_params_hb.view_host()(i,j).a_hb8 = static_cast<KK_FLOAT>(a_hb8[i][j]); k_params_hb.view_host()(j,i).a_hb8 = static_cast<KK_FLOAT>(a_hb8[j][i]);
  k_params_hb.view_host()(i,j).theta_hb8_0 = static_cast<KK_FLOAT>(theta_hb8_0[i][j]); k_params_hb.view_host()(j,i).theta_hb8_0 = static_cast<KK_FLOAT>(theta_hb8_0[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb8_ast = static_cast<KK_FLOAT>(dtheta_hb8_ast[i][j]); k_params_hb.view_host()(j,i).dtheta_hb8_ast = static_cast<KK_FLOAT>(dtheta_hb8_ast[j][i]);
  k_params_hb.view_host()(i,j).b_hb8 = static_cast<KK_FLOAT>(b_hb8[i][j]); k_params_hb.view_host()(j,i).b_hb8 = static_cast<KK_FLOAT>(b_hb8[j][i]);
  k_params_hb.view_host()(i,j).dtheta_hb8_c = static_cast<KK_FLOAT>(dtheta_hb8_c[i][j]); k_params_hb.view_host()(j,i).dtheta_hb8_c = static_cast<KK_FLOAT>(dtheta_hb8_c[j][i]);

  k_params_hb.modify_host();







  // Sync to device
  k_params_hb.template sync<DeviceType>();







  // Register the cutoff of this pair, which includes the displacement of the
  // interaction sites from the COM, with the COM screen of the npair fix,
  // which takes the max over all consuming styles and type pairs.
  if (fix_oxdna_npairKK) fix_oxdna_npairKK->request_screen_cutoff(cutone);

  // "cutone" is "cut_hb_hc[i][j]", sets the master list distance cutoff
  return cutone;

}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   whether this style can be part of the fused hbond + oxdna3/xstk kernel:
   only the screened-pair kernel of the oxDNA3 variant can
------------------------------------------------------------------------- */

template<class DeviceType>
bool PairOxdnaHbondKokkos<DeviceType>::fuse_supported() const
{
  return (execution_space != HostKK) && (oxdnaflag == OXDNA3) && (compute_flag != 0);
}

/* ---------------------------------------------------------------------- */

namespace LAMMPS_NS {
template class PairOxdnaHbondKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdnaHbondKokkos<LMPHostType>;
#endif
}
