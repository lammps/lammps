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
   Contributing authors: Ray Shan (SNL), Stan Moore (SNL),
     Evan Weinberg (NVIDIA)

   Nicholas Curtis (AMD), Leopold Grinberd (AMD), and Gina Sitaraman (AMD):
     - Reduced math overhead: enabled specialized calls (e.g., cbrt for a
         cube root instead of pow) and use power/exponential laws to reduce the
         number of exponentials evaluated, etc.
     - Added blocking to the Torsion and (optionally) BuildLists kernels, to
         reduce thread divergence on GPUs
     - Added preview to BuildLists kernels along with full version
------------------------------------------------------------------------- */

#include "pair_reaxff_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "fix_acks2_reaxff_legacy_kokkos.h"
#include "kokkos.h"
#include "math_const.h"
#include "math_special.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "reaxff_api.h"
#include "reaxff_defs.h"

#include <cmath>

/* ---------------------------------------------------------------------- */

using namespace ReaxFF;

namespace LAMMPS_NS {
using namespace MathConst;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_init_kokkos.hpp"

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairReaxFFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  copymode = 1;

  bocnt = hbcnt = 0;

  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag,vflag,0);

  atomKK->sync(execution_space,datamask_read);
  k_params_sing.template sync<DeviceType>();
  k_params_twbp.template sync<DeviceType>();
  k_params_thbp.template sync<DeviceType>();
  k_params_fbp.template sync<DeviceType>();
  k_params_hbp.template sync<DeviceType>();

  if (eflag_either || vflag_either) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  tag = atomKK->k_tag.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atomKK->nlocal;
  newton_pair = force->newton_pair;

  nn = list->inum;
  NN = atom->nlocal + atom->nghost;

  const int inum = list->inum;
  const int ignum = inum + list->gnum;
  NeighListKokkos<DeviceType>* k_list = static_cast<NeighListKokkos<DeviceType>*>(list);
  d_numneigh = k_list->d_numneigh;
  d_neighbors = k_list->d_neighbors;
  d_ilist = k_list->d_ilist;

  if (acks2_flag) {
    auto ifix = modify->get_fix_by_style("^acks2/reax").front();
    if (ifix->execution_space == Host) {
      auto k_s = ((FixACKS2ReaxFFLegacyKokkos<LMPHostType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    } else {
      auto k_s = ((FixACKS2ReaxFFLegacyKokkos<LMPDeviceType>*) ifix)->get_s();
      k_s.sync<DeviceType>();
      d_s = k_s.view<DeviceType>();
    }
  }

  // allocate duplicated memory
  if (need_dup) {
    dup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(f);
    dup_eatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_eatom);
    dup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_vatom);
  } else {
    ndup_f            = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(f);
    ndup_eatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_eatom);
    ndup_vatom        = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_vatom);
  }

  if (eflag_global) {
    for (int i = 0; i < 14; i++)
      pvector[i] = 0.0;
  }

  EV_FLOAT_REAX ev;
  EV_FLOAT_REAX ev_all;

  // Polarization (self)
  if (neighflag == HALF) {
    if (eflag_global)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALF>>(0,inum),*this,ev);
    else if (eflag_atom)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALF>>(0,inum),*this);
  } else { //if (neighflag == HALFTHREAD) {
    if (eflag_global)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALFTHREAD>>(0,inum),*this,ev);
    else if (eflag_atom)
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputePolar<HALFTHREAD>>(0,inum),*this);
  }
  ev_all += ev;
  pvector[13] = ev.ecoul;

  // LJ + Coulomb
  if (api->control->tabulate) {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALF,0>>(0,inum),*this);
    } else if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALFTHREAD,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTabulatedLJCoulomb<HALFTHREAD,0>>(0,inum),*this);
    }
  } else {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALF,0>>(0,inum),*this);
    } else if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALFTHREAD,1>,Kokkos::LaunchBounds<256,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeLJCoulomb<HALFTHREAD,0>,Kokkos::LaunchBounds<256,1>>(0,inum),*this);
    }
  }
  ev_all += ev;
  pvector[10] = ev.evdwl;
  pvector[11] = ev.ecoul;


  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    allocate_array();
  }

  // allocate duplicated memory
  if (need_dup) {
    dup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_dDeltap_self);
    dup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_total_bo);
  } else {
    ndup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_dDeltap_self);
    ndup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_total_bo);
  }

  // Neighbor lists for bond and hbond

  // try, resize if necessary

  int resize = 1;
  while (resize) {
    resize = 0;

    k_resize_bo.h_view() = 0;
    k_resize_bo.modify<LMPHostType>();
    k_resize_bo.sync<DeviceType>();

    k_resize_hb.h_view() = 0;
    k_resize_hb.modify<LMPHostType>();
    k_resize_hb.sync<DeviceType>();

    // zero
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxZero>(0,nmax),*this);

    if (execution_space == Host) { // CPU
      if (neighflag == HALF)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlocking<HALF>>(0,ignum),*this);
      else if (neighflag == HALFTHREAD)
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlocking<HALFTHREAD>>(0,ignum),*this);
    } else {
      if (list_blocking_flag) {
        if (neighflag == HALF)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlockingPreview<HALF>>(0,ignum),*this);
        else if (neighflag == HALFTHREAD)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfBlockingPreview<HALFTHREAD>>(0,ignum),*this);
      } else {
        if (neighflag == HALF)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfPreview<HALF>>(0,ignum),*this);
        else if (neighflag == HALFTHREAD)
          Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsHalfPreview<HALFTHREAD>>(0,ignum),*this);
      }
    }

    k_resize_bo.modify<DeviceType>();
    k_resize_bo.sync<LMPHostType>();
    int resize_bo = k_resize_bo.h_view();
    if (resize_bo) maxbo = MAX(maxbo+MAX(1,maxbo*0.1),resize_bo);

    k_resize_hb.modify<DeviceType>();
    k_resize_hb.sync<LMPHostType>();
    int resize_hb = k_resize_hb.h_view();
    if (resize_hb) maxhb = MAX(maxhb+MAX(1,maxhb*0.1),resize_hb);

    resize = resize_bo || resize_hb;
    if (resize) {
      allocate_array();
      if (need_dup) {
        dup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_dDeltap_self);
        dup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_total_bo);
      } else {
        ndup_dDeltap_self = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_dDeltap_self);
        ndup_total_bo     = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_total_bo);
      }
    }
  }

  if (execution_space != Host) // GPU
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBuildListsFull>(0,ignum),*this);

  // allocate duplicated memory
  if (need_dup)
    dup_CdDelta = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterDuplicated>(d_CdDelta);
  else
    ndup_CdDelta = Kokkos::Experimental::create_scatter_view<Kokkos::Experimental::ScatterSum, Kokkos::Experimental::ScatterNonDuplicated>(d_CdDelta);

  // reduction over duplicated memory
  if (need_dup)
    Kokkos::Experimental::contribute(d_total_bo, dup_total_bo); // needed in BondOrder1

  // Bond order
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder1>(0,ignum),*this);

  if( api->control->ereaxff_flag ) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder2<1>>(0,ignum),*this);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder3<1>>(0,ignum),*this);
  } else {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder2<0>>(0,ignum),*this);
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxBondOrder3<0>>(0,ignum),*this);
  }

  // Bond energy
  if (neighflag == HALF) {
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALF,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALF,0>>(0,inum),*this);
    ev_all += ev;
    pvector[0] = ev.evdwl;
  } else { //if (neighflag == HALFTHREAD) {
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALFTHREAD,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond1<HALFTHREAD,0>>(0,inum),*this);
    ev_all += ev;
    pvector[0] = ev.evdwl;
  }

  // Multi-body corrections
  if (neighflag == HALF) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti1>(0,inum),*this);
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALF,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALF,0>>(0,inum),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti1>(0,inum),*this);
    if (eflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALFTHREAD,1>>(0,inum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeMulti2<HALFTHREAD,0>>(0,inum),*this);
    ev_all += ev;
  }
  pvector[2] = ev.ereax[0];
  pvector[1] = ev.ereax[1]+ev.ereax[2];
  pvector[3] = 0.0;
  ev_all.evdwl += ev.ereax[0] + ev.ereax[1] + ev.ereax[2];

  int count_angular = 0;
  int count_torsion = 0;

  auto& h_count_angular_torsion = k_count_angular_torsion.h_view;
  h_count_angular_torsion(0) = 0;
  h_count_angular_torsion(1) = 0;
  k_count_angular_torsion.template modify<LMPHostType>();
  k_count_angular_torsion.template sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxCountAngularTorsion<false> >(0,inum),*this);

  k_count_angular_torsion.template modify<DeviceType>();
  k_count_angular_torsion.template sync<LMPHostType>();
  count_angular = h_count_angular_torsion(0);
  count_torsion = h_count_angular_torsion(1);

  if (count_angular > (int)d_angular_pack.extent(0))
    MemKK::realloc_kokkos(d_angular_pack,"reaxff:angular_pack",(int)(count_angular * 1.1),2);
  if (count_torsion > (int)d_torsion_pack.extent(0))
    MemKK::realloc_kokkos(d_torsion_pack,"reaxff:torsion_pack",(int)(count_torsion * 1.1),2);

  // need to zero to re-count
  h_count_angular_torsion(0) = 0;
  h_count_angular_torsion(1) = 0;
  k_count_angular_torsion.template modify<LMPHostType>();
  k_count_angular_torsion.template sync<DeviceType>();

  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxCountAngularTorsion<true>>(0,inum),*this);

  // no need to re-sync count_angular, count_torsion

  // Angular
  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALF,1>>(0,count_angular),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALF,0>>(0,count_angular),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALFTHREAD,1>>(0,count_angular),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeAngularPreprocessed<HALFTHREAD,0>>(0,count_angular),*this);
    ev_all += ev;
  }
  pvector[4] = ev.ereax[3];
  pvector[5] = ev.ereax[4];
  pvector[6] = ev.ereax[5];
  ev_all.evdwl += ev.ereax[3] + ev.ereax[4] + ev.ereax[5];

  // Torsion
  if (neighflag == HALF) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALF,1>>(0,count_torsion),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALF,0>>(0,count_torsion),*this);
    ev_all += ev;
  } else { //if (neighflag == HALFTHREAD) {
    if (evflag)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALFTHREAD,1>>(0,count_torsion),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeTorsionPreprocessed<HALFTHREAD,0>>(0,count_torsion),*this);
    ev_all += ev;
  }
  pvector[8] = ev.ereax[6];
  pvector[9] = ev.ereax[7];
  ev_all.evdwl += ev.ereax[6] + ev.ereax[7];

  // Hydrogen Bond
  if (cut_hbsq > 0.0) {
    if (neighflag == HALF) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALF,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALF,0>>(0,inum),*this);
      ev_all += ev;
    } else { //if (neighflag == HALFTHREAD) {
      if (evflag)
        Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALFTHREAD,1>>(0,inum),*this,ev);
      else
        Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeHydrogen<HALFTHREAD,0>>(0,inum),*this);
      ev_all += ev;
    }
  }
  pvector[7] = ev.ereax[8];
  ev_all.evdwl += ev.ereax[8];

  // reduction over duplicated memory
  if (need_dup) {
    Kokkos::Experimental::contribute(d_dDeltap_self, dup_dDeltap_self); // needed in ComputeBond2
    Kokkos::Experimental::contribute(d_CdDelta, dup_CdDelta); // needed in ComputeBond2
  }

  // Bond force
  if (neighflag == HALF) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxUpdateBond<HALF>>(0,ignum),*this);

    if (vflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALF,1>>(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALF,0>>(0,ignum),*this);
    ev_all += ev;
    pvector[0] += ev.evdwl;
  } else { //if (neighflag == HALFTHREAD) {
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxUpdateBond<HALFTHREAD>>(0,ignum),*this);

    if (vflag_either)
      Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALFTHREAD,1>>(0,ignum),*this,ev);
    else
      Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagPairReaxComputeBond2<HALFTHREAD,0>>(0,ignum),*this);
    ev_all += ev;
    pvector[0] += ev.evdwl;
  }

  // reduction over duplicated memory
  if (need_dup)
    Kokkos::Experimental::contribute(f, dup_f);

  if (eflag_global) {
    eng_vdwl += ev_all.evdwl;
    eng_coul += ev_all.ecoul;
  }

  if (vflag_global) {
    virial[0] += ev_all.v[0];
    virial[1] += ev_all.v[1];
    virial[2] += ev_all.v[2];
    virial[3] += ev_all.v[3];
    virial[4] += ev_all.v[4];
    virial[5] += ev_all.v[5];
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

  if (fixspecies_flag)
    FindBondSpecies();

  copymode = 0;

  // free scatterview memory
  if (need_dup) {
    dup_f            = {};
    dup_eatom        = {};
    dup_vatom        = {};
    dup_dDeltap_self = {};
    dup_total_bo     = {};
    dup_CdDelta      = {};
  } else {
    ndup_f            = {};
    ndup_eatom        = {};
    ndup_vatom        = {};
    ndup_dDeltap_self = {};
    ndup_total_bo     = {};
    ndup_CdDelta      = {};
  }

  d_neighbors = typename AT::t_neighbors_2d();
}

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_ljcoulomb_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_neighbor_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_bondorder_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_angular_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_torsion_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_hydrogen_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_bond_kokkos.hpp"

/* ---------------------------------------------------------------------- */

#include "pair_reaxff_other_kokkos.hpp"

/* ---------------------------------------------------------------------- */

template class PairReaxFFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairReaxFFKokkos<LMPHostType>;
#endif
}
