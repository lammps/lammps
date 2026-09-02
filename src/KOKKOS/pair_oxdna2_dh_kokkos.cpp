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

#include "pair_oxdna2_dh_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "math_const.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna2DhKokkos<DeviceType>::PairOxdna2DhKokkos(LAMMPS *lmp) : PairOxdna2Dh(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA2;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdna2DhKokkos<DeviceType>::~PairOxdna2DhKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2DhKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
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
  qeff = atomKK->k_qeff.template view<DeviceType>();

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

  // loop over neighbors of my atoms for compute functors

  EV_FLOAT ev;

  if (evflag) {
    if (neighflag == HALF) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALF,1,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALF,1,1> >(0,anum),*this,ev);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALF,0,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALF,0,1> >(0,anum),*this,ev);
        }
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALFTHREAD,1,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALFTHREAD,1,1> >(0,anum),*this,ev);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALFTHREAD,0,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALFTHREAD,0,1> >(0,anum),*this,ev);
        }
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,FULL,1,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,FULL,1,1> >(0,anum),*this,ev);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,FULL,0,1> >(0,anum),*this,ev);
        } else {
          Kokkos::parallel_reduce(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,FULL,0,1> >(0,anum),*this,ev);
        }
      }
    }
  } else {
    if (neighflag == HALF) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALF,1,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALF,1,0> >(0,anum),*this);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALF,0,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALF,0,0> >(0,anum),*this);
        }
      }
    } else if (neighflag == HALFTHREAD) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALFTHREAD,1,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALFTHREAD,1,0> >(0,anum),*this);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,HALFTHREAD,0,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,HALFTHREAD,0,0> >(0,anum),*this);
        }
      }
    } else if (neighflag == FULL) {
      if (newton_pair) {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,FULL,1,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,FULL,1,0> >(0,anum),*this);
        }
      } else {
        if (oxdnaflag==OXDNA2) {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXDNA2,FULL,0,0> >(0,anum),*this);
        } else {
          Kokkos::parallel_for(OxdnaRangePolicy<DeviceType, TagPairOxdna2DhCompute<OXRNA2,FULL,0,0> >(0,anum),*this);
        }
      }
    }
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

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2DhKokkos<DeviceType>::operator()(TagPairOxdna2DhCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
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
  KK_FLOAT ra_cs0, ra_cs1, ra_cs2;
  KK_FLOAT rtmp_s0, rtmp_s1, rtmp_s2;

  // vector COM-backbone site a
  if constexpr (OXDNAFLAG==OXDNA2) {
    constexpr KK_FLOAT d_cs_x = -0.34;
    constexpr KK_FLOAT d_cs_y = +0.3408;
    ra_cs0 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,0), d_cs_y*d_ny_xtrct(a,0));
    ra_cs1 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,1), d_cs_y*d_ny_xtrct(a,1));
    ra_cs2 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,2), d_cs_y*d_ny_xtrct(a,2));
  } else {
    // OXRNA2
    constexpr KK_FLOAT d_cs_x = -0.4;
    constexpr KK_FLOAT d_cs_z = +0.2;
    ra_cs0 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,0), d_cs_z*d_nz_xtrct(a,0));
    ra_cs1 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,1), d_cs_z*d_nz_xtrct(a,1));
    ra_cs2 = Kokkos::fma(d_cs_x, d_nx_xtrct(a,2), d_cs_z*d_nz_xtrct(a,2));
  }

  rtmp_s0 = x(a,0) + ra_cs0;
  rtmp_s1 = x(a,1) + ra_cs1;
  rtmp_s2 = x(a,2) + ra_cs2;

  KK_ACC_FLOAT ftmp_a0 = 0.0;
  KK_ACC_FLOAT ftmp_a1 = 0.0;
  KK_ACC_FLOAT ftmp_a2 = 0.0;
  KK_ACC_FLOAT ttmp_a0 = 0.0;
  KK_ACC_FLOAT ttmp_a1 = 0.0;
  KK_ACC_FLOAT ttmp_a2 = 0.0;

  const KK_FLOAT qeff_a = qeff(a);

  const int bnum = d_numneigh(a);

  for (int ib = 0; ib < bnum; ib++) {

    int b = d_neighbors(a,ib);
    const KK_FLOAT factor_lj = special_lj[sbmask(b)];
    if (factor_lj == static_cast<KK_FLOAT>(0.0)) continue;
    b &= NEIGHMASK;
    const int btype = type(b);

    KK_FLOAT rb_cs0, rb_cs1, rb_cs2;
    if constexpr (OXDNAFLAG==OXDNA2) {
      constexpr KK_FLOAT d_cs_x = -0.34;
      constexpr KK_FLOAT d_cs_y = +0.3408;
      rb_cs0 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,0), d_cs_y*d_ny_xtrct(b,0));
      rb_cs1 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,1), d_cs_y*d_ny_xtrct(b,1));
      rb_cs2 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,2), d_cs_y*d_ny_xtrct(b,2));
    } else {
      constexpr KK_FLOAT d_cs_x = -0.4;
      constexpr KK_FLOAT d_cs_z = +0.2;
      rb_cs0 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,0), d_cs_z*d_nz_xtrct(b,0));
      rb_cs1 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,1), d_cs_z*d_nz_xtrct(b,1));
      rb_cs2 = Kokkos::fma(d_cs_x, d_nx_xtrct(b,2), d_cs_z*d_nz_xtrct(b,2));
    }

    const KK_FLOAT delx = rtmp_s0 - x(b,0) - rb_cs0;
    const KK_FLOAT dely = rtmp_s1 - x(b,1) - rb_cs1;
    const KK_FLOAT delz = rtmp_s2 - x(b,2) - rb_cs2;
    const KK_FLOAT rsq = Kokkos::fma(delz, delz, Kokkos::fma(dely, dely, delx * delx));

    if (rsq > d_cutsq_dh_c(atype, btype)) continue; // Note the switch of sign, > vs <=, due to using continue

    const KK_FLOAT qeff_b = qeff(b);
    const KK_FLOAT qeff_dh_pf = d_qeff_dh_pf(atype, btype);
    const KK_FLOAT kappa = d_kappa_dh(atype, btype);
    const KK_FLOAT b_dh = d_b_dh(atype, btype);
    const KK_FLOAT cut_dh_ast = d_cut_dh_ast(atype, btype);
    const KK_FLOAT cut_dh_c = d_cut_dh_c(atype, btype);

    const KK_FLOAT rinv = static_cast<KK_FLOAT>(Kokkos::rsqrt(rsq));
    const KK_FLOAT r = rsq * rinv;
    KK_FLOAT fpair;
    KK_FLOAT evdwl_loc = 0.0;

    if (r <= cut_dh_ast) {
      const KK_FLOAT expterm = expf(-kappa * r);
      fpair = qeff_a * qeff_b * qeff_dh_pf * expterm * (kappa + rinv) * rinv * rinv;

      if constexpr (EVFLAG) {
        evdwl_loc = qeff_a * qeff_b * qeff_dh_pf * expterm * rinv;
      }

    } else {
      const KK_FLOAT delrcut = cut_dh_c - r;
      fpair = 2.0 * qeff_a * qeff_b * b_dh * delrcut * rinv;

      if constexpr (EVFLAG) {
        evdwl_loc = qeff_a * qeff_b * b_dh * delrcut * delrcut; // double negative, so safe to keep delrcut as is
      }
    }

    fpair *= factor_lj;
    if constexpr (EVFLAG) {
      evdwl_loc *= factor_lj;
    }

    const KK_ACC_FLOAT delf0 = static_cast<KK_ACC_FLOAT>(delx * fpair);
    const KK_ACC_FLOAT delf1 = static_cast<KK_ACC_FLOAT>(dely * fpair);
    const KK_ACC_FLOAT delf2 = static_cast<KK_ACC_FLOAT>(delz * fpair);

    // apply force and torque to each of 2 atoms
    ftmp_a0 += delf0;
    ftmp_a1 += delf1;
    ftmp_a2 += delf2;
    const KK_ACC_FLOAT delta0 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(ra_cs1, delf2, -ra_cs2*delf1));
    const KK_ACC_FLOAT delta1 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(ra_cs2, delf0, -ra_cs0*delf2));
    const KK_ACC_FLOAT delta2 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(ra_cs0, delf1, -ra_cs1*delf0));
    ttmp_a0 += delta0;
    ttmp_a1 += delta1;
    ttmp_a2 += delta2;
    if ((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal)) {
      a_f(b,0) -= delf0;
      a_f(b,1) -= delf1;
      a_f(b,2) -= delf2;
      const KK_ACC_FLOAT deltb0 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(rb_cs1, delf2, -rb_cs2*delf1));
      const KK_ACC_FLOAT deltb1 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(rb_cs2, delf0, -rb_cs0*delf2));
      const KK_ACC_FLOAT deltb2 = static_cast<KK_ACC_FLOAT>(Kokkos::fma(rb_cs0, delf1, -rb_cs1*delf0));
      a_torque(b,0) -= deltb0;
      a_torque(b,1) -= deltb1;
      a_torque(b,2) -= deltb2;
    }

    // increment energy and virial
    // NOTE: The virial is calculated on the 'molecular' basis.
    // (see G. Ciccotti and J.P. Ryckaert, Comp. Phys. Rep. 4, 345-392 (1986))

    if constexpr (EVFLAG) {
      ev.evdwl += (((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?1.0:0.5)*evdwl_loc;

      if (vflag_either || eflag_atom) {
        this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR>(ev,a,b,evdwl_loc,\
        delf0,delf1,delf2,x(a,0)-x(b,0), x(a,1)-x(b,1), x(a,2)-x(b,2));
      }
    }
  }
  a_f(a,0) += ftmp_a0;
  a_f(a,1) += ftmp_a1;
  a_f(a,2) += ftmp_a2;
  a_torque(a,0) += ttmp_a0;
  a_torque(a,1) += ttmp_a1;
  a_torque(a,2) += ttmp_a2;
}

template<class DeviceType>
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdna2DhKokkos<DeviceType>::operator()(TagPairOxdna2DhCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, \
  const int &ia) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>\
  (TagPairOxdna2DhCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>(),ia,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2DhKokkos<DeviceType>::allocate()
{
  PairOxdna2Dh::allocate();

  int n = atom->ntypes;

  memoryKK->create_kokkos(k_qeff_dh_pf,n+1,n+1,"PairOxdna2Dh:qeff_dh_pf");
  memoryKK->create_kokkos(k_kappa_dh,n+1,n+1,"PairOxdna2Dh:kappa_dh");
  memoryKK->create_kokkos(k_b_dh,n+1,n+1,"PairOxdna2Dh:b_dh");
  memoryKK->create_kokkos(k_cut_dh_ast,n+1,n+1,"PairOxdna2Dh:cut_dh_ast");
  memoryKK->create_kokkos(k_cutsq_dh_ast,n+1,n+1,"PairOxdna2Dh:cutsq_dh_ast");
  memoryKK->create_kokkos(k_cut_dh_c,n+1,n+1,"PairOxdna2Dh:cut_dh_c");
  memoryKK->create_kokkos(k_cutsq_dh_c,n+1,n+1,"PairOxdna2Dh:cutsq_dh_c");

  d_qeff_dh_pf = k_qeff_dh_pf.template view<DeviceType>();
  d_kappa_dh = k_kappa_dh.template view<DeviceType>();
  d_b_dh = k_b_dh.template view<DeviceType>();
  d_cut_dh_ast = k_cut_dh_ast.template view<DeviceType>();
  d_cutsq_dh_ast = k_cutsq_dh_ast.template view<DeviceType>();
  d_cut_dh_c = k_cut_dh_c.template view<DeviceType>();
  d_cutsq_dh_c = k_cutsq_dh_c.template view<DeviceType>();

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2DhKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdna2DhKokkos<DeviceType>::init_style()
{
  neighbor->add_request(this);
  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();

  fix_oxdna_lrfKK = nullptr;
  Kokkos::fence("before OXDNA/LRF/kk lookup");
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdna2DhKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdna2Dh::init_one(i,j);

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  k_qeff_dh_pf.view_host()(i,j) = qeff_dh_pf[i][j]; k_qeff_dh_pf.view_host()(j,i) = qeff_dh_pf[j][i];
  k_kappa_dh.view_host()(i,j) = kappa_dh[i][j]; k_kappa_dh.view_host()(j,i) = kappa_dh[j][i];
  k_b_dh.view_host()(i,j) = b_dh[i][j]; k_b_dh.view_host()(j,i) = b_dh[j][i];
  k_cut_dh_ast.view_host()(i,j) = cut_dh_ast[i][j]; k_cut_dh_ast.view_host()(j,i) = cut_dh_ast[j][i];
  k_cutsq_dh_ast.view_host()(i,j) = cutsq_dh_ast[i][j]; k_cutsq_dh_ast.view_host()(j,i) = cutsq_dh_ast[j][i];
  k_cut_dh_c.view_host()(i,j) = cut_dh_c[i][j]; k_cut_dh_c.view_host()(j,i) = cut_dh_c[j][i];
  k_cutsq_dh_c.view_host()(i,j) = cutsq_dh_c[i][j]; k_cutsq_dh_c.view_host()(j,i) = cutsq_dh_c[j][i];

  k_qeff_dh_pf.template modify<LMPHostType>();
  k_kappa_dh.template modify<LMPHostType>();
  k_b_dh.template modify<LMPHostType>();
  k_cut_dh_ast.template modify<LMPHostType>();
  k_cutsq_dh_ast.template modify<LMPHostType>();
  k_cut_dh_c.template modify<LMPHostType>();
  k_cutsq_dh_c.template modify<LMPHostType>();

  // Sync to device
  k_qeff_dh_pf.template sync<DeviceType>();
  k_kappa_dh.template sync<DeviceType>();
  k_b_dh.template sync<DeviceType>();
  k_cut_dh_ast.template sync<DeviceType>();
  k_cutsq_dh_ast.template sync<DeviceType>();
  k_cut_dh_c.template sync<DeviceType>();
  k_cutsq_dh_c.template sync<DeviceType>();

  // "cutone" is "cut_dh_c[i][j]", sets the master list distance cutoff
  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR>
KOKKOS_INLINE_FUNCTION
void PairOxdna2DhKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
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
int PairOxdna2DhKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

namespace LAMMPS_NS {
template class PairOxdna2DhKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdna2DhKokkos<LMPHostType>;
#endif
}
