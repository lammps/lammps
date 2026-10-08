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

#include "pair_oxdna_excv_kokkos.h"

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
#include "fix_oxdna_prime_neighs_kokkos.h"
#include "mf_oxdna_kokkos.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaExcvKokkos<DeviceType>::PairOxdnaExcvKokkos(LAMMPS *lmp) : PairOxdnaExcv(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA;
  fix_oxdna_lrfKK = nullptr;
  fix_oxdna_npairKK = nullptr;
  fix_oxdna_prime_neighsKK = nullptr;
  last_prime_neighs_atom_nbuild = -1;
  params2_uniform = 0;
  params2_dirty = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaExcvKokkos<DeviceType>::~PairOxdnaExcvKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }

  if (fix_oxdna_lrfKK) modify->delete_fix(fix_oxdna_lrfKK->id);

  // also remove the other internal helper fixes, so they do not keep requesting
  // neighbor lists after the oxDNA styles are gone. Styles that still need them
  // create them again in their init_style().

  if (modify->get_fix_by_id("npair_kk")) modify->delete_fix("npair_kk");
  if (modify->get_fix_by_id("prime_neighs_kk")) modify->delete_fix("prime_neighs_kk");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  xn_type = fix_oxdna_lrfKK->packed_type();
  xn = fix_oxdna_lrfKK->packed();
  f = atomKK->k_f.template view<DeviceType>();
  torque = atomKK->k_torque.template view<DeviceType>();
  type = atomKK->k_type.template view<DeviceType>();

  nlocal = atom->nlocal;
  newton_pair = force->newton_pair;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];

  tag = atomKK->k_tag.template view<DeviceType>();
  id5p = atomKK->k_id5p.template view<DeviceType>();
  id3p = atomKK->k_id3p.template view<DeviceType>();

  map_style = atom->map_style;
  if (map_style == Atom::MAP_ARRAY) {
    k_map_array = atomKK->k_map_array;
    k_map_array.template sync<DeviceType>();
  } else if (map_style == Atom::MAP_HASH) {
    k_map_hash = atomKK->k_map_hash;
    k_map_hash.template sync<DeviceType>();
  }

  // types of the 3'/5' neighbors of all atoms, for the base-base terms of
  // bonded pairs; recomputed after every reneighboring

  if (neighbor->nbuild != last_prime_neighs_atom_nbuild) {
    fix_oxdna_prime_neighsKK->compute_prime_neighs_atom();
    last_prime_neighs_atom_nbuild = neighbor->nbuild;
  }
  d_prime_neighs_atom = fix_oxdna_prime_neighsKK->d_prime_neighs_atom;

  // get the neighbor list and neighbors used in operator()

  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_neighbors = k_list->d_neighbors;
  anum = list->inum;
  d_alist = k_list->d_ilist;
  d_numneigh = k_list->d_numneigh;

  // split the neighbors of an atom over several threads on GPUs when there
  // are too few atoms to fill a quarter of the device
  // (only with atomic updates of atom a, i.e. a half list with HALFTHREAD)
  nsplit = 1;
  if ((execution_space != HostKK) && (neighflag == HALFTHREAD)) {
    const int target = DeviceType().concurrency() / 4;
    while ((nsplit < 4) && (anum * nsplit < target)) nsplit *= 2;
  }

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

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->packed_nx();
  d_ny_xtrct = fix_oxdna_lrfKK->packed_ny();
  d_nz_xtrct = fix_oxdna_lrfKK->packed_nz();

  // check whether the non-tetramer coefficients are the same for all type pairs

  if (params2_dirty) {
    const auto h_params2 = k_params2_excv.view_host();
    const int n = atom->ntypes;
    params2_uniform = 1;
    for (int i = 1; i <= n; i++)
      for (int j = 1; j <= n; j++)
        if (std::memcmp(&h_params2(i,j), &h_params2(1,1), sizeof(ParamsOxdnaExcv2)) != 0)
          params2_uniform = 0;
    params2_uni = h_params2(1,1);
    params2_dirty = 0;
  }

  // loop over neighbors of my atoms for compute functors

  copymode = 1;

  EV_FLOAT ev;

  auto run_compute_by_oxdnaflag = [&](auto neighflag_tag, auto newtonpair_tag, auto evflag_tag) {
    constexpr int NEIGHFLAG = decltype(neighflag_tag)::value;
    constexpr int NEWTON_PAIR = decltype(newtonpair_tag)::value;
    constexpr int EVFLAG = decltype(evflag_tag)::value;

    auto run_compute = [&](auto model_tag) {
      constexpr int MODEL = decltype(model_tag)::value;
      if constexpr (EVFLAG) {
        Kokkos::parallel_reduce(
          OxdnaRangePolicy<DeviceType, TagPairOxdnaExcvCompute<MODEL,NEIGHFLAG,NEWTON_PAIR,EVFLAG> >(0,anum*nsplit),
          *this,ev);
      } else {
        Kokkos::parallel_for(
          OxdnaRangePolicy<DeviceType, TagPairOxdnaExcvCompute<MODEL,NEIGHFLAG,NEWTON_PAIR,EVFLAG> >(0,anum*nsplit),
          *this);
      }
    };

    if (oxdnaflag == OXDNA) {
      run_compute(std::integral_constant<int,OXDNA>{});
    } else if (oxdnaflag == OXDNA2) {
      run_compute(std::integral_constant<int,OXDNA2>{});
    } else if (oxdnaflag == OXRNA2) {
      run_compute(std::integral_constant<int,OXRNA2>{});
    } else if (oxdnaflag == OXDNA3) {
      run_compute(std::integral_constant<int,OXDNA3>{});
    } else {
      error->all(FLERR, "Unknown OXDNA model flag in pair oxdna/excv/kk");
    }
  };

  const int dispatch_neigh =
      (neighflag == HALF) ? 0 :
      (neighflag == HALFTHREAD) ? 1 :
      (neighflag == FULL) ? 2 : -1;

  if (dispatch_neigh < 0) {
    error->all(FLERR, "Unsupported neighbor flag in pair oxdna/excv/kk");
  }

  const int dispatch_key = (evflag ? 8 : 0) | (newton_pair ? 4 : 0) | dispatch_neigh;
  switch (dispatch_key) {
    case 0: run_compute_by_oxdnaflag(std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 1: run_compute_by_oxdnaflag(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 2: run_compute_by_oxdnaflag(std::integral_constant<int,FULL>{},       std::integral_constant<int,0>{}, std::integral_constant<int,0>{}); break;
    case 4: run_compute_by_oxdnaflag(std::integral_constant<int,HALF>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 5: run_compute_by_oxdnaflag(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 6: run_compute_by_oxdnaflag(std::integral_constant<int,FULL>{},       std::integral_constant<int,1>{}, std::integral_constant<int,0>{}); break;
    case 8: run_compute_by_oxdnaflag(std::integral_constant<int,HALF>{},       std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 9: run_compute_by_oxdnaflag(std::integral_constant<int,HALFTHREAD>{}, std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 10: run_compute_by_oxdnaflag(std::integral_constant<int,FULL>{},      std::integral_constant<int,0>{}, std::integral_constant<int,1>{}); break;
    case 12: run_compute_by_oxdnaflag(std::integral_constant<int,HALF>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    case 13: run_compute_by_oxdnaflag(std::integral_constant<int,HALFTHREAD>{},std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    case 14: run_compute_by_oxdnaflag(std::integral_constant<int,FULL>{},      std::integral_constant<int,1>{}, std::integral_constant<int,1>{}); break;
    default: error->all(FLERR, "Internal dispatch error in pair oxdna/excv/kk");
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

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairOxdnaExcvKokkos<DeviceType>::map_tag(const tagint &itag) const
{
  int mapped = -1;
  if (itag == -1) return mapped;
  if (map_style == Atom::MAP_ARRAY) {
    const auto map_array = k_map_array.view<DeviceType>();
    if ((itag >= 0) && (itag < static_cast<tagint>(map_array.extent(0)))) mapped = map_array(itag);
  } else if (map_style == Atom::MAP_HASH) {
    mapped = AtomKokkos::map_find_hash_kokkos<DeviceType>(itag, k_map_hash);
  }
  return mapped;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaExcvKokkos<DeviceType>::operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia, EV_FLOAT &ev) const
{
  // f and torque array are duplicated for OpenMP, atomic for GPU, and neither for Serial

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<AtomicDup_v<NEIGHFLAG,DeviceType>>();

  // with nsplit > 1, each atom's neighbors are shared by nsplit consecutive
  // work items (more threads for small systems); a-side updates are atomic
  // anyway, so this only changes the summation order
  const int a = d_alist(ia / nsplit);
  const int ib0 = ia % nsplit;
  const int atype = type(a);
  // vectors COM-backbone site in lab frame
  KK_FLOAT ra_cbk[3], rb_cbk[3];
  KK_FLOAT ra_cbs[3], rb_cbs[3];
  KK_FLOAT rtmp_bk[3], rtmp_bs[3];

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];    // force, torque increment
  KK_ACC_FLOAT evdwl, fpair;                   // energy, force
  KK_FLOAT delr_bkbk[3],rsq_bkbk,delr_bkbs[3],rsq_bkbs;
  KK_FLOAT delr_bs[3],rsq_bs,delr_bsbs[3],rsq_bsbs;

  KK_ACC_FLOAT ftmp[3],ttmp[3];  // temporary atom-a force and torque to reduce excessive dup/atomic updates.

  // vector COM - backbone and base sites a
  if constexpr (OXDNAFLAG==OXDNA) {
    constexpr KK_FLOAT dx_cbk_oxdna1=static_cast<KK_FLOAT>(-0.4);
    ra_cbk[0] = dx_cbk_oxdna1*d_nx_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxdna1*d_nx_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxdna1*d_nx_xtrct(a,2);
    constexpr KK_FLOAT dx_cbs_oxdna1=static_cast<KK_FLOAT>(+0.4);
    ra_cbs[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_cbs[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_cbs[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
  } else if constexpr (OXDNAFLAG==OXDNA2) {
    constexpr KK_FLOAT dx_cbk_oxdna2 = static_cast<KK_FLOAT>(-0.34);
    constexpr KK_FLOAT dy_cbk_oxdna2 = static_cast<KK_FLOAT>(+0.3408);
    ra_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(a,0) + dy_cbk_oxdna2*d_ny_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(a,1) + dy_cbk_oxdna2*d_ny_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(a,2) + dy_cbk_oxdna2*d_ny_xtrct(a,2);
    constexpr KK_FLOAT dx_cbs_oxdna1 = static_cast<KK_FLOAT>(+0.4);  // same base sites as OXDNA1
    ra_cbs[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_cbs[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_cbs[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
  } else if constexpr (OXDNAFLAG==OXRNA2) {
    constexpr KK_FLOAT dx_cbk_oxrna2 = static_cast<KK_FLOAT>(-0.4);
    constexpr KK_FLOAT dz_cbk_oxrna2 = static_cast<KK_FLOAT>(+0.2);
    ra_cbk[0] = dx_cbk_oxrna2*d_nx_xtrct(a,0) + dz_cbk_oxrna2*d_nz_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxrna2*d_nx_xtrct(a,1) + dz_cbk_oxrna2*d_nz_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxrna2*d_nx_xtrct(a,2) + dz_cbk_oxrna2*d_nz_xtrct(a,2);
    constexpr KK_FLOAT dx_cbs_oxdna1 = static_cast<KK_FLOAT>(+0.4);  // same base sites as OXDNA1
    ra_cbs[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_cbs[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_cbs[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
  } else if constexpr (OXDNAFLAG==OXDNA3) {
    // Uses the same backbone sites as OXDNA2...
    constexpr KK_FLOAT dx_cbk_oxdna2 = static_cast<KK_FLOAT>(-0.34);
    constexpr KK_FLOAT dy_cbk_oxdna2 = static_cast<KK_FLOAT>(+0.3408);
    ra_cbk[0] = dx_cbk_oxdna2*d_nx_xtrct(a,0) + dy_cbk_oxdna2*d_ny_xtrct(a,0);
    ra_cbk[1] = dx_cbk_oxdna2*d_nx_xtrct(a,1) + dy_cbk_oxdna2*d_ny_xtrct(a,1);
    ra_cbk[2] = dx_cbk_oxdna2*d_nx_xtrct(a,2) + dy_cbk_oxdna2*d_ny_xtrct(a,2);
    // ...but different base sites than OXDNA2.
    constexpr KK_FLOAT dx_cbs_pur_oxdna3 = static_cast<KK_FLOAT>(+0.43);
    constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = static_cast<KK_FLOAT>(+0.37);
    int nucl_acid = (atype%4);
    if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
      ra_cbs[0] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,0);
      ra_cbs[1] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,1);
      ra_cbs[2] = dx_cbs_pyr_oxdna3*d_nx_xtrct(a,2);
    } else {  // purine (A or G)
      ra_cbs[0] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,0);
      ra_cbs[1] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,1);
      ra_cbs[2] = dx_cbs_pur_oxdna3*d_nx_xtrct(a,2);
    }
  }

  rtmp_bk[0] = x(a,0)+ra_cbk[0];
  rtmp_bk[1] = x(a,1)+ra_cbk[1];
  rtmp_bk[2] = x(a,2)+ra_cbk[2];
  rtmp_bs[0] = x(a,0)+ra_cbs[0];
  rtmp_bs[1] = x(a,1)+ra_cbs[1];
  rtmp_bs[2] = x(a,2)+ra_cbs[2];

  const int bnum = d_numneigh(a);

  // topology of atom a, used to find its bonded (base-base) partners
  const tagint tag_a = tag(a);
  const tagint id3p_a = id3p(a);
  const tagint id5p_a = id5p(a);

  ftmp[0] = 0.0;
  ftmp[1] = 0.0;
  ftmp[2] = 0.0;
  ttmp[0] = 0.0;
  ttmp[1] = 0.0;
  ttmp[2] = 0.0;

  for (int ib = ib0; ib < bnum; ib += nsplit) {

    int b = d_neighbors(a,ib);
    const KK_FLOAT factor_lj = static_cast<KK_FLOAT>(special_lj[sbmask(b)]);
    // the 3'/5' neighbors of a are 1-2 special neighbors, whose special bits
    // this list keeps; only those need the topology loads below
    const bool bonded12 = (sbmask(b) == 1);
    b &= NEIGHMASK;
    // the packed record of b with 16-byte loads: position and type first,
    // the frame vectors only for neighbors within the center-of-mass cutoff
    OxdnaRow rowb;
    oxdna_load_row<4>(xn, b, rowb);
    const int btype = static_cast<int>(rowb.v[3]);

    // no sites interact beyond this center-of-mass distance
    const KK_FLOAT dx_com = x(a,0) - rowb.v[0];
    const KK_FLOAT dy_com = x(a,1) - rowb.v[1];
    const KK_FLOAT dz_com = x(a,2) - rowb.v[2];
    if (dx_com*dx_com + dy_com*dy_com + dz_com*dz_com > d_cutsq_com(atype,btype)) continue;
    oxdna_load_row_rest<4,16>(xn, b, rowb);

    const ParamsOxdnaExcv2 p2 = params2(atype,btype);

    // b-side force and torque, summed over the site pairs of this neighbor
    // and applied with one atomic update per component
    KK_ACC_FLOAT fb[3] = {0.0, 0.0, 0.0};
    KK_ACC_FLOAT tb[3] = {0.0, 0.0, 0.0};

    // vector COM - backbone and base sites b
    if constexpr (OXDNAFLAG==OXDNA) {
      constexpr KK_FLOAT dx_cbk_oxdna1=static_cast<KK_FLOAT>(-0.4);
      rb_cbk[0] = dx_cbk_oxdna1*rowb.v[4];
      rb_cbk[1] = dx_cbk_oxdna1*rowb.v[5];
      rb_cbk[2] = dx_cbk_oxdna1*rowb.v[6];
      constexpr KK_FLOAT dx_cbs_oxdna1=static_cast<KK_FLOAT>(+0.4);
      rb_cbs[0] = dx_cbs_oxdna1*rowb.v[4];
      rb_cbs[1] = dx_cbs_oxdna1*rowb.v[5];
      rb_cbs[2] = dx_cbs_oxdna1*rowb.v[6];
    } else if constexpr (OXDNAFLAG==OXDNA2) {
      constexpr KK_FLOAT dx_cbk_oxdna2 = static_cast<KK_FLOAT>(-0.34);
      constexpr KK_FLOAT dy_cbk_oxdna2 = static_cast<KK_FLOAT>(+0.3408);
      rb_cbk[0] = dx_cbk_oxdna2*rowb.v[4] + dy_cbk_oxdna2*rowb.v[7];
      rb_cbk[1] = dx_cbk_oxdna2*rowb.v[5] + dy_cbk_oxdna2*rowb.v[8];
      rb_cbk[2] = dx_cbk_oxdna2*rowb.v[6] + dy_cbk_oxdna2*rowb.v[9];
      constexpr KK_FLOAT dx_cbs_oxdna1 = static_cast<KK_FLOAT>(+0.4);  // same base sites as OXDNA1
      rb_cbs[0] = dx_cbs_oxdna1*rowb.v[4];
      rb_cbs[1] = dx_cbs_oxdna1*rowb.v[5];
      rb_cbs[2] = dx_cbs_oxdna1*rowb.v[6];
    } else if constexpr (OXDNAFLAG==OXRNA2) {
      constexpr KK_FLOAT dx_cbk_oxrna2 = static_cast<KK_FLOAT>(-0.4);
      constexpr KK_FLOAT dz_cbk_oxrna2 = static_cast<KK_FLOAT>(+0.2);
      rb_cbk[0] = dx_cbk_oxrna2*rowb.v[4] + dz_cbk_oxrna2*rowb.v[10];
      rb_cbk[1] = dx_cbk_oxrna2*rowb.v[5] + dz_cbk_oxrna2*rowb.v[11];
      rb_cbk[2] = dx_cbk_oxrna2*rowb.v[6] + dz_cbk_oxrna2*rowb.v[12];
      constexpr KK_FLOAT dx_cbs_oxdna1 = static_cast<KK_FLOAT>(+0.4);  // same base sites as OXDNA1
      rb_cbs[0] = dx_cbs_oxdna1*rowb.v[4];
      rb_cbs[1] = dx_cbs_oxdna1*rowb.v[5];
      rb_cbs[2] = dx_cbs_oxdna1*rowb.v[6];
    } else if constexpr (OXDNAFLAG==OXDNA3) {
      // Uses the same backbone sites as OXDNA2...
      constexpr KK_FLOAT dx_cbk_oxdna2 = static_cast<KK_FLOAT>(-0.34);
      constexpr KK_FLOAT dy_cbk_oxdna2 = static_cast<KK_FLOAT>(+0.3408);
      rb_cbk[0] = dx_cbk_oxdna2*rowb.v[4] + dy_cbk_oxdna2*rowb.v[7];
      rb_cbk[1] = dx_cbk_oxdna2*rowb.v[5] + dy_cbk_oxdna2*rowb.v[8];
      rb_cbk[2] = dx_cbk_oxdna2*rowb.v[6] + dy_cbk_oxdna2*rowb.v[9];
      // ...but different base sites than OXDNA2.
      constexpr KK_FLOAT dx_cbs_pur_oxdna3 = static_cast<KK_FLOAT>(+0.43);
      constexpr KK_FLOAT dx_cbs_pyr_oxdna3 = static_cast<KK_FLOAT>(+0.37);
      int nucl_acid = (btype%4);
      if (nucl_acid==0 || nucl_acid==2) {  // pyrimdine (C or T)
        rb_cbs[0] = dx_cbs_pyr_oxdna3*rowb.v[4];
        rb_cbs[1] = dx_cbs_pyr_oxdna3*rowb.v[5];
        rb_cbs[2] = dx_cbs_pyr_oxdna3*rowb.v[6];
      } else {  // purine (A or G)
        rb_cbs[0] = dx_cbs_pur_oxdna3*rowb.v[4];
        rb_cbs[1] = dx_cbs_pur_oxdna3*rowb.v[5];
        rb_cbs[2] = dx_cbs_pur_oxdna3*rowb.v[6];
      }
    }

    // vector backbone site b to a
    delr_bkbk[0] = rtmp_bk[0] - (rowb.v[0]+rb_cbk[0]);
    delr_bkbk[1] = rtmp_bk[1] - (rowb.v[1]+rb_cbk[1]);
    delr_bkbk[2] = rtmp_bk[2] - (rowb.v[2]+rb_cbk[2]);
    rsq_bkbk = delr_bkbk[0]*delr_bkbk[0] + delr_bkbk[1]*delr_bkbk[1] + delr_bkbk[2]*delr_bkbk[2];
    // vector base site b to backbone site a
    delr_bkbs[0] = rtmp_bk[0] - (rowb.v[0]+rb_cbs[0]);
    delr_bkbs[1] = rtmp_bk[1] - (rowb.v[1]+rb_cbs[1]);
    delr_bkbs[2] = rtmp_bk[2] - (rowb.v[2]+rb_cbs[2]);
    rsq_bkbs = delr_bkbs[0]*delr_bkbs[0] + delr_bkbs[1]*delr_bkbs[1] + delr_bkbs[2]*delr_bkbs[2];
    // vector backbone site b to base site a
    delr_bs[0] = rtmp_bs[0] - (rowb.v[0]+rb_cbk[0]);
    delr_bs[1] = rtmp_bs[1] - (rowb.v[1]+rb_cbk[1]);
    delr_bs[2] = rtmp_bs[2] - (rowb.v[2]+rb_cbk[2]);
    rsq_bs = delr_bs[0]*delr_bs[0] + delr_bs[1]*delr_bs[1] + delr_bs[2]*delr_bs[2];
    // vector base site b to a
    delr_bsbs[0] = rtmp_bs[0] - (rowb.v[0]+rb_cbs[0]);
    delr_bsbs[1] = rtmp_bs[1] - (rowb.v[1]+rb_cbs[1]);
    delr_bsbs[2] = rtmp_bs[2] - (rowb.v[2]+rb_cbs[2]);
    rsq_bsbs = delr_bsbs[0]*delr_bsbs[0] + delr_bsbs[1]*delr_bsbs[1] + delr_bsbs[2]*delr_bsbs[2];

    // excluded volume interactions:

    // backbone-backbone
    if (rsq_bkbk < p2.cutsq_bkbk_c) {
      // F3 modulation factor, force and energy calculation
      evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bkbk,p2.cutsq_bkbk_ast,p2.cut_bkbk_c,p2.lj1_bkbk,
                        p2.lj2_bkbk,p2.epsilon_bkbk,p2.b_bkbk,fpair));
      // knock out nearest-neighbor interaction between ss
      fpair *= static_cast<KK_ACC_FLOAT>(factor_lj);
      evdwl *= static_cast<KK_ACC_FLOAT>(factor_lj);
      // force and torque increment calculation
      delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbk[0]);
      delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbk[1]);
      delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbk[2]);
      delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[1];
      delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[2];
      delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[0];
      ftmp[0] += delf[0];
      ftmp[1] += delf[1];
      ftmp[2] += delf[2];
      ttmp[0] += delta[0];
      ttmp[1] += delta[1];
      ttmp[2] += delta[2];
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
        fb[0] -= delf[0];
        fb[1] -= delf[1];
        fb[2] -= delf[2];
        deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[1];
        deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[2];
        deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[0];
        tb[0] -= deltb[0];
        tb[1] -= deltb[1];
        tb[2] -= deltb[2];
      }
      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
        }

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
          delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
        }
      }
    }

    // backbone-base
    if (rsq_bkbs < p2.cutsq_bkbs_c) {
      // F3 modulation factor, force and energy calculation
      evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bkbs,p2.cutsq_bkbs_ast,p2.cut_bkbs_c,p2.lj1_bkbs,
                        p2.lj2_bkbs,p2.epsilon_bkbs,p2.b_bkbs,fpair));
      // force and torque increment calculation
      delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbs[0]);
      delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbs[1]);
      delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bkbs[2]);
      delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[1];
      delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[2];
      delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[0];
      ftmp[0] += delf[0];
      ftmp[1] += delf[1];
      ftmp[2] += delf[2];
      ttmp[0] += delta[0];
      ttmp[1] += delta[1];
      ttmp[2] += delta[2];
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
        fb[0] -= delf[0];
        fb[1] -= delf[1];
        fb[2] -= delf[2];
        deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[1];
        deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[2];
        deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[0];
        tb[0] -= deltb[0];
        tb[1] -= deltb[1];
        tb[2] -= deltb[2];
      }
      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
        }

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
          delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
        }
      }
    }

    // base-backbone
    if (rsq_bs < p2.cutsq_bkbs_c) {
      // F3 modulation factor, force and energy calculation
      evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bs,p2.cutsq_bkbs_ast,p2.cut_bkbs_c,p2.lj1_bkbs,
                        p2.lj2_bkbs,p2.epsilon_bkbs,p2.b_bkbs,fpair));
      // force and torque increment calculation
      delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bs[0]);
      delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bs[1]);
      delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bs[2]);
      delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[1];
      delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[2];
      delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[0];
      ftmp[0] += delf[0];
      ftmp[1] += delf[1];
      ftmp[2] += delf[2];
      ttmp[0] += delta[0];
      ttmp[1] += delta[1];
      ttmp[2] += delta[2];
      if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
        fb[0] -= delf[0];
        fb[1] -= delf[1];
        fb[2] -= delf[2];
        deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[1];
        deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[2];
        deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[0];
        tb[0] -= deltb[0];
        tb[1] -= deltb[1];
        tb[2] -= deltb[2];
      }
      if (EVFLAG) {
        if (eflag) {
          ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
        }

        if (vflag_either || eflag_atom) {
          this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
          delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
        }
      }
    }

    // base-base
    if (bonded12 && (tag_a == id3p(b)) && (tag(b) == id5p_a)) {
      // types of the 3' neighbor of a and the 5' neighbor of b (0 for a strand end)
      const int _3ptype = d_prime_neighs_atom(a,2);
      const int _5ptype = d_prime_neighs_atom(b,3);
      if (rsq_bsbs < d_params4_excv(_3ptype,atype,btype,_5ptype).cut4sq_bsbs_c) {
        // F3 modulation factor, force and energy calculation
        evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bsbs,d_params4_excv(_3ptype,atype,btype,_5ptype).cut4sq_bsbs_ast,d_params4_excv(_3ptype,atype,btype,_5ptype).cut4_bsbs_c,
                          d_params4_excv(_3ptype,atype,btype,_5ptype).lj14_bsbs,d_params4_excv(_3ptype,atype,btype,_5ptype).lj24_bsbs,
                          p2.epsilon_bsbs,d_params4_excv(_3ptype,atype,btype,_5ptype).b4_bsbs,fpair));
        // force and torque increment calculation
        delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[0]);
        delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[1]);
        delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[2]);
        delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[1];
        delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[2];
        delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[0];
        ftmp[0] += delf[0];
        ftmp[1] += delf[1];
        ftmp[2] += delf[2];
        ttmp[0] += delta[0];
        ttmp[1] += delta[1];
        ttmp[2] += delta[2];
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
          fb[0] -= delf[0];
          fb[1] -= delf[1];
          fb[2] -= delf[2];
          deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[1];
          deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[2];
          deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[0];
          tb[0] -= deltb[0];
          tb[1] -= deltb[1];
          tb[2] -= deltb[2];
        }
        if (EVFLAG) {
          if (eflag) {
            ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
          }

          if (vflag_either || eflag_atom) {
            this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
            delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
          }
        }
      }
    } else if (bonded12 && (tag_a == id5p(b)) && (tag(b) == id3p_a)) {
      // types of the 3' neighbor of b and the 5' neighbor of a (0 for a strand end)
      const int _3ptype = d_prime_neighs_atom(b,2);
      const int _5ptype = d_prime_neighs_atom(a,3);
      if (rsq_bsbs < d_params4_excv(_3ptype,btype,atype,_5ptype).cut4sq_bsbs_c) {
        // F3 modulation factor, force and energy calculation
        evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bsbs,d_params4_excv(_3ptype,btype,atype,_5ptype).cut4sq_bsbs_ast,d_params4_excv(_3ptype,btype,atype,_5ptype).cut4_bsbs_c,
                          d_params4_excv(_3ptype,btype,atype,_5ptype).lj14_bsbs,d_params4_excv(_3ptype,btype,atype,_5ptype).lj24_bsbs,
                          params2(btype,atype).epsilon_bsbs,d_params4_excv(_3ptype,btype,atype,_5ptype).b4_bsbs,fpair));
        // force and torque increment calculation
        delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[0]);
        delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[1]);
        delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[2]);
        delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[1];
        delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[2];
        delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[0];
        ftmp[0] += delf[0];
        ftmp[1] += delf[1];
        ftmp[2] += delf[2];
        ttmp[0] += delta[0];
        ttmp[1] += delta[1];
        ttmp[2] += delta[2];
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
          fb[0] -= delf[0];
          fb[1] -= delf[1];
          fb[2] -= delf[2];
          deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[1];
          deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[2];
          deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[0];
          tb[0] -= deltb[0];
          tb[1] -= deltb[1];
          tb[2] -= deltb[2];
        }
        if (EVFLAG) {
          if (eflag) {
            ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
          }

          if (vflag_either || eflag_atom) {
            this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
            delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
          }
        }
      }
    } else {
      if (rsq_bsbs < p2.cutsq_bsbs_c) {
        // F3 modulation factor, force and energy calculation
        evdwl = static_cast<KK_ACC_FLOAT>(F3_KK(rsq_bsbs,p2.cutsq_bsbs_ast,p2.cut_bsbs_c,p2.lj1_bsbs,
                          p2.lj2_bsbs,p2.epsilon_bsbs,p2.b_bsbs,fpair));
        // force and torque increment calculation
        delf[0] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[0]);
        delf[1] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[1]);
        delf[2] = fpair * static_cast<KK_ACC_FLOAT>(delr_bsbs[2]);
        delta[0] = static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[1];
        delta[1] = static_cast<KK_ACC_FLOAT>(ra_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[2];
        delta[2] = static_cast<KK_ACC_FLOAT>(ra_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbs[1])*delf[0];
        ftmp[0] += delf[0];
        ftmp[1] += delf[1];
        ftmp[2] += delf[2];
        ttmp[0] += delta[0];
        ttmp[1] += delta[1];
        ttmp[2] += delta[2];
        if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
          fb[0] -= delf[0];
          fb[1] -= delf[1];
          fb[2] -= delf[2];
          deltb[0] = static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[1];
          deltb[1] = static_cast<KK_ACC_FLOAT>(rb_cbs[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[2];
          deltb[2] = static_cast<KK_ACC_FLOAT>(rb_cbs[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbs[1])*delf[0];
          tb[0] -= deltb[0];
          tb[1] -= deltb[1];
          tb[2] -= deltb[2];
        }
        if (EVFLAG) {
          if (eflag) {
            ev.evdwl += static_cast<KK_ACC_FLOAT>(((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl;
          }

          if (vflag_either || eflag_atom) {
            this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,static_cast<KK_FLOAT>(evdwl),\
            delf[0],delf[1],delf[2],x(a,0)-rowb.v[0], x(a,1)-rowb.v[1], x(a,2)-rowb.v[2]);
          }
        }
      }
    }
    // end excluded volume interaction
    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
      if ((fb[0] != 0.0) || (fb[1] != 0.0) || (fb[2] != 0.0) ||
          (tb[0] != 0.0) || (tb[1] != 0.0) || (tb[2] != 0.0)) {
        a_f(b,0) += fb[0];
        a_f(b,1) += fb[1];
        a_f(b,2) += fb[2];
        a_torque(b,0) += tb[0];
        a_torque(b,1) += tb[1];
        a_torque(b,2) += tb[2];
      }
    }
  }
  a_f(a,0) += ftmp[0];
  a_f(a,1) += ftmp[1];
  a_f(a,2) += ftmp[2];
  a_torque(a,0) += ttmp[0];
  a_torque(a,1) += ttmp[1];
  a_torque(a,2) += ttmp[2];
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaExcvKokkos<DeviceType>::operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ia,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::allocate()
{
  PairOxdnaExcv::allocate();

  int n = atom->ntypes;

  k_params2_excv = decltype(k_params2_excv)("PairOxdnaExcvKokkos:params2_excv", n+1, n+1);
  k_params4_excv = decltype(k_params4_excv)("PairOxdnaExcvKokkos:params4_excv", n+1, n+1, n+1, n+1);
  k_cutsq_com = decltype(k_cutsq_com)("PairOxdnaExcvKokkos:cutsq_com", n+1, n+1);
  d_params2_excv = k_params2_excv.template view<DeviceType>();
  d_params4_excv = k_params4_excv.template view<DeviceType>();
  d_cutsq_com = k_cutsq_com.template view<DeviceType>();
  params2_dirty = 1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR, "The oxDNA and oxRNA pair styles do not take any arguments");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::init_style()
{
  // the internal helper fixes are always created for the default KOKKOS variant,
  // so /kk/host styles cannot work with them when LAMMPS is compiled for a GPU

  if (std::is_same_v<DeviceType, LMPHostType> && !std::is_same_v<DeviceType, LMPDeviceType>)
    error->all(FLERR, "The /kk/host variants of the CG-DNA styles are not supported "
               "when LAMMPS is compiled for a GPU");

  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  // ensure fix OXDNA/LRF/kk is added for backward-compatability
  if (!fix_oxdna_lrfKK) {
    fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(modify->add_fix("lrf_kk all OXDNA/LRF/kk"));
  }
  // ensure fix OXDNA/NPAIR/kk exists; reuse an existing one, since adding a fix
  // with the same ID would delete the instance other styles refer to
  auto npair_fixes = modify->get_fix_by_style("^OXDNA/NPAIR/kk");
  if (npair_fixes.size() == 0)
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(modify->add_fix("npair_kk all OXDNA/NPAIR/kk"));
  else
    fix_oxdna_npairKK = dynamic_cast<FixOxdnaNpairKokkos<DeviceType> *>(npair_fixes[0]);
  if (!fix_oxdna_npairKK) error->all(FLERR, "Fix OXDNA/NPAIR/kk not found");

  auto prime_fixes = modify->get_fix_by_style("^OXDNA/PRIME_NEIGHS/kk");
  if (prime_fixes.size() == 0) {
    fix_oxdna_prime_neighsKK =
      dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(modify->add_fix("prime_neighs_kk all OXDNA/PRIME_NEIGHS/kk"));
  } else {
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(prime_fixes[0]);
  }

  if (!fix_oxdna_prime_neighsKK) error->all(FLERR, "Fix OXDNA/PRIME_NEIGHS/kk not found");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdnaExcvKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdnaExcv::init_one(i,j);

  // All non-tetramer Kokkos views are set here within ::init_one, and
  // the tetramer Kokkos views are set within ::coeff

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  k_params2_excv.view_host()(i,j).epsilon_bkbk = static_cast<KK_FLOAT>(epsilon_bkbk[i][j]);
  k_params2_excv.view_host()(j,i).epsilon_bkbk = static_cast<KK_FLOAT>(epsilon_bkbk[j][i]);
  k_params2_excv.view_host()(i,j).sigma_bkbk = static_cast<KK_FLOAT>(sigma_bkbk[i][j]);
  k_params2_excv.view_host()(j,i).sigma_bkbk = static_cast<KK_FLOAT>(sigma_bkbk[j][i]);
  k_params2_excv.view_host()(i,j).cut_bkbk_ast = static_cast<KK_FLOAT>(cut_bkbk_ast[i][j]);
  k_params2_excv.view_host()(j,i).cut_bkbk_ast = static_cast<KK_FLOAT>(cut_bkbk_ast[j][i]);
  k_params2_excv.view_host()(i,j).b_bkbk = static_cast<KK_FLOAT>(b_bkbk[i][j]);
  k_params2_excv.view_host()(j,i).b_bkbk = static_cast<KK_FLOAT>(b_bkbk[j][i]);
  k_params2_excv.view_host()(i,j).cut_bkbk_c = static_cast<KK_FLOAT>(cut_bkbk_c[i][j]);
  k_params2_excv.view_host()(j,i).cut_bkbk_c = static_cast<KK_FLOAT>(cut_bkbk_c[j][i]);
  k_params2_excv.view_host()(i,j).lj1_bkbk = static_cast<KK_FLOAT>(lj1_bkbk[i][j]);
  k_params2_excv.view_host()(j,i).lj1_bkbk = static_cast<KK_FLOAT>(lj1_bkbk[j][i]);
  k_params2_excv.view_host()(i,j).lj2_bkbk = static_cast<KK_FLOAT>(lj2_bkbk[i][j]);
  k_params2_excv.view_host()(j,i).lj2_bkbk = static_cast<KK_FLOAT>(lj2_bkbk[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bkbk_ast = static_cast<KK_FLOAT>(cutsq_bkbk_ast[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bkbk_ast = static_cast<KK_FLOAT>(cutsq_bkbk_ast[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bkbk_c = static_cast<KK_FLOAT>(cutsq_bkbk_c[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bkbk_c = static_cast<KK_FLOAT>(cutsq_bkbk_c[j][i]);

  k_params2_excv.view_host()(i,j).epsilon_bkbs = static_cast<KK_FLOAT>(epsilon_bkbs[i][j]);
  k_params2_excv.view_host()(j,i).epsilon_bkbs = static_cast<KK_FLOAT>(epsilon_bkbs[j][i]);
  k_params2_excv.view_host()(i,j).sigma_bkbs = static_cast<KK_FLOAT>(sigma_bkbs[i][j]);
  k_params2_excv.view_host()(j,i).sigma_bkbs = static_cast<KK_FLOAT>(sigma_bkbs[j][i]);
  k_params2_excv.view_host()(i,j).cut_bkbs_ast = static_cast<KK_FLOAT>(cut_bkbs_ast[i][j]);
  k_params2_excv.view_host()(j,i).cut_bkbs_ast = static_cast<KK_FLOAT>(cut_bkbs_ast[j][i]);
  k_params2_excv.view_host()(i,j).b_bkbs = static_cast<KK_FLOAT>(b_bkbs[i][j]);
  k_params2_excv.view_host()(j,i).b_bkbs = static_cast<KK_FLOAT>(b_bkbs[j][i]);
  k_params2_excv.view_host()(i,j).cut_bkbs_c = static_cast<KK_FLOAT>(cut_bkbs_c[i][j]);
  k_params2_excv.view_host()(j,i).cut_bkbs_c = static_cast<KK_FLOAT>(cut_bkbs_c[j][i]);
  k_params2_excv.view_host()(i,j).lj1_bkbs = static_cast<KK_FLOAT>(lj1_bkbs[i][j]);
  k_params2_excv.view_host()(j,i).lj1_bkbs = static_cast<KK_FLOAT>(lj1_bkbs[j][i]);
  k_params2_excv.view_host()(i,j).lj2_bkbs = static_cast<KK_FLOAT>(lj2_bkbs[i][j]);
  k_params2_excv.view_host()(j,i).lj2_bkbs = static_cast<KK_FLOAT>(lj2_bkbs[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bkbs_ast = static_cast<KK_FLOAT>(cutsq_bkbs_ast[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bkbs_ast = static_cast<KK_FLOAT>(cutsq_bkbs_ast[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bkbs_c = static_cast<KK_FLOAT>(cutsq_bkbs_c[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bkbs_c = static_cast<KK_FLOAT>(cutsq_bkbs_c[j][i]);

  k_params2_excv.view_host()(i,j).epsilon_bsbs = static_cast<KK_FLOAT>(epsilon_bsbs[i][j]);
  k_params2_excv.view_host()(j,i).epsilon_bsbs = static_cast<KK_FLOAT>(epsilon_bsbs[j][i]);
  k_params2_excv.view_host()(i,j).sigma_bsbs = static_cast<KK_FLOAT>(sigma_bsbs[i][j]);
  k_params2_excv.view_host()(j,i).sigma_bsbs = static_cast<KK_FLOAT>(sigma_bsbs[j][i]);
  k_params2_excv.view_host()(i,j).cut_bsbs_ast = static_cast<KK_FLOAT>(cut_bsbs_ast[i][j]);
  k_params2_excv.view_host()(j,i).cut_bsbs_ast = static_cast<KK_FLOAT>(cut_bsbs_ast[j][i]);
  k_params2_excv.view_host()(i,j).b_bsbs = static_cast<KK_FLOAT>(b_bsbs[i][j]);
  k_params2_excv.view_host()(j,i).b_bsbs = static_cast<KK_FLOAT>(b_bsbs[j][i]);
  k_params2_excv.view_host()(i,j).cut_bsbs_c = static_cast<KK_FLOAT>(cut_bsbs_c[i][j]);
  k_params2_excv.view_host()(j,i).cut_bsbs_c = static_cast<KK_FLOAT>(cut_bsbs_c[j][i]);
  k_params2_excv.view_host()(i,j).lj1_bsbs = static_cast<KK_FLOAT>(lj1_bsbs[i][j]);
  k_params2_excv.view_host()(j,i).lj1_bsbs = static_cast<KK_FLOAT>(lj1_bsbs[j][i]);
  k_params2_excv.view_host()(i,j).lj2_bsbs = static_cast<KK_FLOAT>(lj2_bsbs[i][j]);
  k_params2_excv.view_host()(j,i).lj2_bsbs = static_cast<KK_FLOAT>(lj2_bsbs[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bsbs_ast = static_cast<KK_FLOAT>(cutsq_bsbs_ast[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bsbs_ast = static_cast<KK_FLOAT>(cutsq_bsbs_ast[j][i]);
  k_params2_excv.view_host()(i,j).cutsq_bsbs_c = static_cast<KK_FLOAT>(cutsq_bsbs_c[i][j]);
  k_params2_excv.view_host()(j,i).cutsq_bsbs_c = static_cast<KK_FLOAT>(cutsq_bsbs_c[j][i]);

  k_params2_excv.modify_host();
  params2_dirty = 1;

  // the margin of 0.01 (in units of the cutoff) keeps rounding of the site
  // positions from ever skipping an interacting pair
  const double cutsq_com = (cutone + 0.01) * (cutone + 0.01);
  k_cutsq_com.view_host()(i,j) = static_cast<KK_FLOAT>(cutsq_com);
  k_cutsq_com.view_host()(j,i) = static_cast<KK_FLOAT>(cutsq_com);
  k_cutsq_com.modify_host();

  // Sync to device
  k_params2_excv.template sync<DeviceType>();
  k_cutsq_com.template sync<DeviceType>();

  // "cutone" is "cut_bkbk_c[i][j]", sets the master list distance cutoff
  return cutone;

}

/* ----------------------------------------------------------------------
   Helper function to set the tetramer Kokkos views within ::coeff
   Is used within child stk/kk classes too.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::coeff_set_tetramers_kokkos(int narg, char **arg)
{
  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          k_params4_excv.view_host()(i,j,k,l).sigma4_bsbs = static_cast<KK_FLOAT>(sigma4_bsbs[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).cut4_bsbs_ast = static_cast<KK_FLOAT>(cut4_bsbs_ast[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).cut4sq_bsbs_ast = static_cast<KK_FLOAT>(cut4sq_bsbs_ast[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).lj14_bsbs = static_cast<KK_FLOAT>(lj14_bsbs[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).lj24_bsbs = static_cast<KK_FLOAT>(lj24_bsbs[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).b4_bsbs = static_cast<KK_FLOAT>(b4_bsbs[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).cut4_bsbs_c = static_cast<KK_FLOAT>(cut4_bsbs_c[i][j][k][l]);
          k_params4_excv.view_host()(i,j,k,l).cut4sq_bsbs_c = static_cast<KK_FLOAT>(cut4sq_bsbs_c[i][j][k][l]);
        }
      }
    }
  }

  k_params4_excv.modify_host();

  // Sync to device
  k_params4_excv.template sync<DeviceType>();
}

/* ----------------------------------------------------------------------
   The tetramer Kokkos views are set here within ::coeff, and the
   non-tetramer Kokkos views are set within ::init_one
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaExcvKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairOxdnaExcv::coeff(narg,arg);

  coeff_set_tetramers_kokkos(narg,arg);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxdnaExcvKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
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
    const KK_ACC_FLOAT v0 = static_cast<KK_ACC_FLOAT>(delx)*fx;
    const KK_ACC_FLOAT v1 = static_cast<KK_ACC_FLOAT>(dely)*fy;
    const KK_ACC_FLOAT v2 = static_cast<KK_ACC_FLOAT>(delz)*fz;
    const KK_ACC_FLOAT v3 = static_cast<KK_ACC_FLOAT>(delx)*fy;
    const KK_ACC_FLOAT v4 = static_cast<KK_ACC_FLOAT>(delx)*fz;
    const KK_ACC_FLOAT v5 = static_cast<KK_ACC_FLOAT>(dely)*fz;

    if (vflag_global) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5) * v0;
          ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5) * v1;
          ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5) * v2;
          ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5) * v3;
          ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5) * v4;
          ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5) * v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5) * v0;
        ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5) * v1;
        ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5) * v2;
        ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5) * v3;
        ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5) * v4;
        ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5) * v5;
        }
      } else {
        ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5) * v0;
        ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5) * v1;
        ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5) * v2;
        ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5) * v3;
        ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5) * v4;
        ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5) * v5;
      }
    }

    if (vflag_atom) {
      if (NEIGHFLAG!=FULL) {
        if (NEWTON_PAIR || i < nlocal) {
          a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5) * v0;
          a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5) * v1;
          a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5) * v2;
          a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5) * v3;
          a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5) * v4;
          a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5) * v5;
        }
        if (NEWTON_PAIR || j < nlocal) {
        a_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5) * v0;
        a_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5) * v1;
        a_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5) * v2;
        a_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5) * v3;
        a_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5) * v4;
        a_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5) * v5;
        }
      } else {
        a_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5) * v0;
        a_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5) * v1;
        a_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5) * v2;
        a_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5) * v3;
        a_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5) * v4;
        a_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5) * v5;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
int PairOxdnaExcvKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairOxdnaExcvKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdnaExcvKokkos<LMPHostType>;
#endif
}
