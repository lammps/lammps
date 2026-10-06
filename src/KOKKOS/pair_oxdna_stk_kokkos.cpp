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

#include "pair_oxdna_stk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "modify.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_prime_neighs_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaStkKokkos<DeviceType>::PairOxdnaStkKokkos(LAMMPS *lmp) : PairOxdnaStk(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  oxdnaflag = EnabledOXDNAFlag::OXDNA;
  fix_oxdna_prime_neighsKK = nullptr;
  last_prime_neighs_bond_nbuild = -1;
  tetramer_uniform = 0;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxdnaStkKokkos<DeviceType>::~PairOxdnaStkKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

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
  xn = fix_oxdna_lrfKK->packed();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;
  nbondlist = neighborKK->nbondlist;

  // Keep bond-context precompute aligned with the current neighbor-list epoch.
  if (last_prime_neighs_bond_nbuild != neighbor->nbuild) {
    fix_oxdna_prime_neighsKK->compute_prime_neighs_bond(d_prime_neighs_bond_own);
    last_prime_neighs_bond_nbuild = neighbor->nbuild;
  }

  d_prime_neighs_bond = d_prime_neighs_bond_own;

  int need_dup = lmp->kokkos->need_dup<DeviceType>();

  copymode = 1;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->packed_nx();
  d_ny_xtrct = fix_oxdna_lrfKK->packed_ny();
  d_nz_xtrct = fix_oxdna_lrfKK->packed_nz();

  // loop over neighbors of my atoms for compute functors

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      if (oxdnaflag==OXDNA) {
        Kokkos::parallel_reduce(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA,1,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag==OXDNA3) {
        Kokkos::parallel_reduce(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA3,1,1> >(0,nbondlist),*this,ev);
      }
    } else {
      if (oxdnaflag==OXDNA) {
        Kokkos::parallel_reduce(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA,0,1> >(0,nbondlist),*this,ev);
      } else if (oxdnaflag==OXDNA3) {
        Kokkos::parallel_reduce(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA3,0,1> >(0,nbondlist),*this,ev);
      }
    }
  } else {
    if (newton_bond) {
      if (oxdnaflag==OXDNA) {
        Kokkos::parallel_for(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA,1,0> >(0,nbondlist),*this);
      } else if (oxdnaflag==OXDNA3) {
        Kokkos::parallel_for(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA3,1,0> >(0,nbondlist),*this);
      }
    } else {
      if (oxdnaflag==OXDNA) {
        Kokkos::parallel_for(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA,0,0> >(0,nbondlist),*this);
      } else if (oxdnaflag==OXDNA3) {
        Kokkos::parallel_for(OxdnaBondRangePolicy<DeviceType, TagPairOxdnaStkCompute<OXDNA3,0,0> >(0,nbondlist),*this);
      }
    }
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
    dup_eatom    = decltype(dup_eatom)();
    dup_vatom    = decltype(dup_vatom)();
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaStkKokkos<DeviceType>::operator()(TagPairOxdnaStkCompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, \
  const int &in, EV_FLOAT &ev) const
{
  // The f and torque arrays are updated with Kokkos::atomic_add() on the views
  // directly (an atomic-trait view copy made here would live in local memory)

  // Use precomputed bond and prime neighbors.
  // NOTE: already in correct order from precompute, so directionality test: a -> b is 3' -> 5' is already satisfied
  int a = d_prime_neighs_bond(in,0);
  int b = d_prime_neighs_bond(in,1);
  // packed records of the atoms with 16-byte loads
  OxdnaRow rowa;
  oxdna_load_row<16>(xn, a, rowa);
  OxdnaRow rowb;
  oxdna_load_row<16>(xn, b, rowb);

  int a3ptype,atype,btype,b5ptype;

  KK_FLOAT ra_cstk[3], rb_cstk[3];           // vectors COM-stacking sites in lab frame
  KK_FLOAT ra_cbk[3], rb_cbk[3];             // vectors COM-backbone sites in lab frame

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];    // force, torque increment
  KK_ACC_FLOAT fsum[3], tsuma[3], tsumb[3];    // total force and torques of this bond
  KK_ACC_FLOAT evdwl,finc,tpair;
  KK_FLOAT delr_bkbk[3],delr_bkbk_norm[3],rsq_bkbk,r_bkbk,rinv_bkbk;
  KK_FLOAT delr_stkstk[3],delr_stkstk_norm[3],rsq_stkstk,r_stkstk,rinv_stkstk;
  KK_FLOAT theta4,t4dir[3],cost4;
  KK_FLOAT theta5p,t5pdir[3],cost5p;
  KK_FLOAT theta6p,t6pdir[3],cost6p;
  KK_FLOAT cosphi1,cosphi2,cosphi1dir[3],cosphi2dir[3];

  KK_FLOAT f1,f4t4,f4t5,f4t6,f5c1,f5c2;
  KK_FLOAT df1,df4t4,df4t5,df4t6,df5c1,df5c2;

  // vector COM [a/b] - stacking site [a/b]
  if constexpr (OXDNAFLAG==OXDNA) {
    // Used for oxDNA[1] and oxDNA2, but not oxDNA3
    constexpr KK_FLOAT dx_cstk_oxdna1 = static_cast<KK_FLOAT>(+0.34);
    ra_cstk[0] = dx_cstk_oxdna1 * rowa.v[4];
    ra_cstk[1] = dx_cstk_oxdna1 * rowa.v[5];
    ra_cstk[2] = dx_cstk_oxdna1 * rowa.v[6];
    rb_cstk[0] = dx_cstk_oxdna1 * rowb.v[4];
    rb_cstk[1] = dx_cstk_oxdna1 * rowb.v[5];
    rb_cstk[2] = dx_cstk_oxdna1 * rowb.v[6];
  } else if constexpr (OXDNAFLAG==OXDNA3) {
    constexpr KK_FLOAT dx_cstk_oxdna3 = static_cast<KK_FLOAT>(+0.37);
    ra_cstk[0] = dx_cstk_oxdna3 * rowa.v[4];
    ra_cstk[1] = dx_cstk_oxdna3 * rowa.v[5];
    ra_cstk[2] = dx_cstk_oxdna3 * rowa.v[6];
    rb_cstk[0] = dx_cstk_oxdna3 * rowb.v[4];
    rb_cstk[1] = dx_cstk_oxdna3 * rowb.v[5];
    rb_cstk[2] = dx_cstk_oxdna3 * rowb.v[6];
  }

  // vector stacking site a to b
  delr_stkstk[0] = rowb.v[0] + rb_cstk[0] - rowa.v[0] - ra_cstk[0];
  delr_stkstk[1] = rowb.v[1] + rb_cstk[1] - rowa.v[1] - ra_cstk[1];
  delr_stkstk[2] = rowb.v[2] + rb_cstk[2] - rowa.v[2] - ra_cstk[2];

  // determine tetramer types
  // Our prime_neighs_bond ordering (a,b,id3p[a],id5p[b]) from precompute
  // is assigned such that we preserve the vanilla oxDNA convention of:
  // 3'neighbor a - a - b - 5'neighbor b
  // throughout the rest of compute.
  // If none of the tetramer-indexed coefficients depend on the 3'/5' context
  // types, index them with 0,0 and skip the context lookups (bit-identical).
  if (tetramer_uniform) {
    a3ptype = 0;
    b5ptype = 0;
  } else {
    const int id3p_local = d_prime_neighs_bond(in,2);
    a3ptype = (id3p_local != -1) ? type(id3p_local) : 0;
    const int id5p_local = d_prime_neighs_bond(in,3);
    b5ptype = (id5p_local != -1) ? type(id5p_local) : 0;
  }

  atype = type(a);
  btype = type(b);

  rsq_stkstk = Kokkos::fma(delr_stkstk[0], delr_stkstk[0], Kokkos::fma(delr_stkstk[1], delr_stkstk[1], delr_stkstk[2]*delr_stkstk[2]));
  r_stkstk = Kokkos::sqrt(rsq_stkstk);
  rinv_stkstk = static_cast<KK_FLOAT>(1.0)/r_stkstk;

  delr_stkstk_norm[0] = delr_stkstk[0] * rinv_stkstk;
  delr_stkstk_norm[1] = delr_stkstk[1] * rinv_stkstk;
  delr_stkstk_norm[2] = delr_stkstk[2] * rinv_stkstk;

  // vector COM [a/b] - backbone site [a/b]
  // All oxDNA variants use the same COM-backbone site offset, so we can use a single constexpr here.
  constexpr KK_FLOAT dx_cbk_oxdna = static_cast<KK_FLOAT>(-0.4);
  ra_cbk[0] = dx_cbk_oxdna * rowa.v[4];
  ra_cbk[1] = dx_cbk_oxdna * rowa.v[5];
  ra_cbk[2] = dx_cbk_oxdna * rowa.v[6];
  rb_cbk[0] = dx_cbk_oxdna * rowb.v[4];
  rb_cbk[1] = dx_cbk_oxdna * rowb.v[5];
  rb_cbk[2] = dx_cbk_oxdna * rowb.v[6];

  // vector backbone site a to b
  delr_bkbk[0] = rowb.v[0] + rb_cbk[0] - rowa.v[0] - ra_cbk[0];
  delr_bkbk[1] = rowb.v[1] + rb_cbk[1] - rowa.v[1] - ra_cbk[1];
  delr_bkbk[2] = rowb.v[2] + rb_cbk[2] - rowa.v[2] - ra_cbk[2];

  rsq_bkbk = Kokkos::fma(delr_bkbk[0], delr_bkbk[0], Kokkos::fma(delr_bkbk[1], delr_bkbk[1], delr_bkbk[2]*delr_bkbk[2]));
  r_bkbk = Kokkos::sqrt(rsq_bkbk);
  rinv_bkbk = static_cast<KK_FLOAT>(1.0)/r_bkbk;

  delr_bkbk_norm[0] = delr_bkbk[0] * rinv_bkbk;
  delr_bkbk_norm[1] = delr_bkbk[1] * rinv_bkbk;
  delr_bkbk_norm[2] = delr_bkbk[2] * rinv_bkbk;

  // beginning of modulation factors

  // f1 = f1 modulation factor
  f1 = F1_KK(r_stkstk, d_params2_st(atype, btype).epsilon_st, d_params2_st(atype, btype).a_st, d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_0,
          d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_lc, d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_hc,
          d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_lo, d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_hi,
          d_params2_st(atype, btype).b_st_lo, d_params2_st(atype, btype).b_st_hi, d_params4_st(a3ptype,atype,btype,b5ptype).shift_st);

  // start early rejection criterium
  if (f1 == static_cast<KK_FLOAT>(0.0)) return;

  // theta4 angle and correction
  cost4 = rowb.v[10] * rowa.v[10] +
          rowb.v[11] * rowa.v[11] +
          rowb.v[12] * rowa.v[12];
  if (cost4 > static_cast<KK_FLOAT>(1.0)) cost4 = static_cast<KK_FLOAT>(1.0);
  if (cost4 < static_cast<KK_FLOAT>(-1.0)) cost4 = static_cast<KK_FLOAT>(-1.0);
  // sin(theta) from the cross product and theta from atan2 stay accurate
  // near 0 and pi, and avoid the slow path of sin(acos(x))
  const KK_FLOAT nz_a[3] = {rowa.v[10], rowa.v[11], rowa.v[12]};
  const KK_FLOAT nz_b[3] = {rowb.v[10], rowb.v[11], rowb.v[12]};
  const KK_FLOAT sin4 = cross_norm(nz_b, nz_a);
  theta4 = Kokkos::atan2(sin4, cost4);
  // f4t4 = f4 modulation factor
  f4t4 = F4_KK(theta4, d_params4_st(a3ptype,atype,btype,b5ptype).a_st4, d_params2_st(atype, btype).theta_st4_0,
               d_params4_st(a3ptype,atype,btype,b5ptype).dtheta_st4_ast, d_params4_st(a3ptype,atype,btype,b5ptype).b_st4,
               d_params4_st(a3ptype,atype,btype,b5ptype).dtheta_st4_c);

  // early rejection criterium
  if (f4t4 == static_cast<KK_FLOAT>(0.0)) return;

  // theta5 angle and correction
  cost5p = rowb.v[10] * delr_stkstk_norm[0] +
           rowb.v[11] * delr_stkstk_norm[1] +
           rowb.v[12] * delr_stkstk_norm[2];
  if (cost5p > static_cast<KK_FLOAT>(1.0)) cost5p = static_cast<KK_FLOAT>(1.0);
  if (cost5p < static_cast<KK_FLOAT>(-1.0)) cost5p = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT sin5p = cross_norm(nz_b, delr_stkstk_norm);
  theta5p = Kokkos::atan2(sin5p, cost5p);
  // f4t5 = f4 modulation factor
  f4t5 = F4_KK(theta5p, d_params2_st(atype, btype).a_st5, d_params2_st(atype, btype).theta_st5_0,
               d_params2_st(atype, btype).dtheta_st5_ast, d_params2_st(atype, btype).b_st5, d_params2_st(atype, btype).dtheta_st5_c);

  // early rejection criterium
  if (f4t5 == static_cast<KK_FLOAT>(0.0)) return;

  // theta6 angle and correction
  cost6p = delr_stkstk_norm[0] * rowa.v[10] +
           delr_stkstk_norm[1] * rowa.v[11] +
           delr_stkstk_norm[2] * rowa.v[12];
  if (cost6p > static_cast<KK_FLOAT>(1.0)) cost6p = static_cast<KK_FLOAT>(1.0);
  if (cost6p < static_cast<KK_FLOAT>(-1.0)) cost6p = static_cast<KK_FLOAT>(-1.0);
  const KK_FLOAT sin6p = cross_norm(delr_stkstk_norm, nz_a);
  theta6p = Kokkos::atan2(sin6p, cost6p);
  // cosphi1 and cosphi2 angles
  cosphi1 = delr_bkbk_norm[0] * rowb.v[7] +
            delr_bkbk_norm[1] * rowb.v[8] +
            delr_bkbk_norm[2] * rowb.v[9];
  cosphi2 = delr_bkbk_norm[0] * rowa.v[7] +
            delr_bkbk_norm[1] * rowa.v[8] +
            delr_bkbk_norm[2] * rowa.v[9];
  if (cosphi1 > static_cast<KK_FLOAT>(1.0)) cosphi1 = static_cast<KK_FLOAT>(1.0);
  if (cosphi1 < static_cast<KK_FLOAT>(-1.0)) cosphi1 = static_cast<KK_FLOAT>(-1.0);
  if (cosphi2 > static_cast<KK_FLOAT>(1.0)) cosphi2 = static_cast<KK_FLOAT>(1.0);
  if (cosphi2 < static_cast<KK_FLOAT>(-1.0)) cosphi2 = static_cast<KK_FLOAT>(-1.0);
  // f4t6 = f4 modulation factor
  f4t6 = F4_KK(theta6p, d_params2_st(atype, btype).a_st6, d_params2_st(atype, btype).theta_st6_0,
               d_params2_st(atype, btype).dtheta_st6_ast, d_params2_st(atype, btype).b_st6, d_params2_st(atype, btype).dtheta_st6_c);
  // f5c1 = f5 modulation factor
  f5c1 = F5_KK(-cosphi1, d_params2_st(atype, btype).a_st1, -d_params2_st(atype, btype).cosphi_st1_ast,
               d_params2_st(atype, btype).b_st1, -d_params2_st(atype, btype).cosphi_st1_c);
  // f5c2 = f5 modulation factor
  f5c2 = F5_KK(-cosphi2, d_params2_st(atype, btype).a_st2, -d_params2_st(atype, btype).cosphi_st2_ast,
               d_params2_st(atype, btype).b_st2, -d_params2_st(atype, btype).cosphi_st2_c);

  evdwl = static_cast<KK_ACC_FLOAT>(f1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2);

  // early rejection criterium
  if (evdwl == static_cast<KK_ACC_FLOAT>(0.0)) return;

  // df1 = derivative of f1 modulation factor
  df1 = DF1_KK(r_stkstk, d_params2_st(atype, btype).epsilon_st, d_params2_st(atype, btype).a_st,
      d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_0,
      d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_lc, d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_hc,
      d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_lo, d_params4_st(a3ptype,atype,btype,b5ptype).cut_st_hi,
      d_params2_st(atype, btype).b_st_lo, d_params2_st(atype, btype).b_st_hi);
  // df4t4 = derivative of f4 modulation factor
  df4t4 = DF4_KK(theta4, d_params4_st(a3ptype,atype,btype,b5ptype).a_st4, d_params2_st(atype, btype).theta_st4_0,
      d_params4_st(a3ptype,atype,btype,b5ptype).dtheta_st4_ast, d_params4_st(a3ptype,atype,btype,b5ptype).b_st4,
      d_params4_st(a3ptype,atype,btype,b5ptype).dtheta_st4_c)/sin4;
  // df4t5 = derivative of f4 modulation factor
  df4t5 = DF4_KK(theta5p, d_params2_st(atype, btype).a_st5, d_params2_st(atype, btype).theta_st5_0, d_params2_st(atype, btype).dtheta_st5_ast,
      d_params2_st(atype, btype).b_st5, d_params2_st(atype, btype).dtheta_st5_c)/sin5p;
  // df4t6 = derivative of f4 modulation factor
  df4t6 = DF4_KK(theta6p, d_params2_st(atype, btype).a_st6, d_params2_st(atype, btype).theta_st6_0, d_params2_st(atype, btype).dtheta_st6_ast,
      d_params2_st(atype, btype).b_st6, d_params2_st(atype, btype).dtheta_st6_c)/sin6p;
  // df5c1 = derivative of f5 modulation factor
  df5c1 = DF5_KK(-cosphi1, d_params2_st(atype, btype).a_st1, -d_params2_st(atype, btype).cosphi_st1_ast,
      d_params2_st(atype, btype).b_st1, -d_params2_st(atype, btype).cosphi_st1_c);
  // df5c2 = derivative of f5 modulation factor
  df5c2 = DF5_KK(-cosphi2, d_params2_st(atype, btype).a_st2, -d_params2_st(atype, btype).cosphi_st2_ast,
      d_params2_st(atype, btype).b_st2, -d_params2_st(atype, btype).cosphi_st2_c);

  // force, torque and virial contribution for forces between stacking sites
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
  finc = static_cast<KK_ACC_FLOAT>(-df1 * f4t4 * f4t5 * f4t6 * f5c1 * f5c2);

  delf[0] += static_cast<KK_ACC_FLOAT>(delr_stkstk[0]) * finc;
  delf[1] += static_cast<KK_ACC_FLOAT>(delr_stkstk[1]) * finc;
  delf[2] += static_cast<KK_ACC_FLOAT>(delr_stkstk[2]) * finc;

  // theta5p force
  if (theta5p != static_cast<KK_FLOAT>(0.0)) {
    finc = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2 * rinv_stkstk);

    delf[0] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[0]*cost5p - rowb.v[10]) * finc;
    delf[1] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[1]*cost5p - rowb.v[11]) * finc;
    delf[2] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[2]*cost5p - rowb.v[12]) * finc;
  }

  // theta6p force
  if (theta6p != static_cast<KK_FLOAT>(0.0)) {
    finc = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2 * rinv_stkstk);

    delf[0] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[0]*cost6p - rowa.v[10]) * finc;
    delf[1] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[1]*cost6p - rowa.v[11]) * finc;
    delf[2] += static_cast<KK_ACC_FLOAT>(delr_stkstk_norm[2]*cost6p - rowa.v[12]) * finc;
  }

  // accumulate the force and torques of this bond; they are applied with a
  // single atomic update per atom and component at the end
  fsum[0] = delf[0];
  fsum[1] = delf[1];
  fsum[2] = delf[2];
  tsuma[0] = static_cast<KK_ACC_FLOAT>(ra_cstk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cstk[2])*delf[1];
  tsuma[1] = static_cast<KK_ACC_FLOAT>(ra_cstk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cstk[0])*delf[2];
  tsuma[2] = static_cast<KK_ACC_FLOAT>(ra_cstk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cstk[1])*delf[0];
  tsumb[0] = static_cast<KK_ACC_FLOAT>(rb_cstk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cstk[2])*delf[1];
  tsumb[1] = static_cast<KK_ACC_FLOAT>(rb_cstk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cstk[0])*delf[2];
  tsumb[2] = static_cast<KK_ACC_FLOAT>(rb_cstk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cstk[1])*delf[0];

  if (EVFLAG) { ev_tally_xyz(ev, a, b, nlocal, NEWTON_BOND, static_cast<KK_FLOAT>(evdwl), delf[0], delf[1], delf[2], \
    rowb.v[0]-rowa.v[0], rowb.v[1]-rowa.v[1], rowb.v[2]-rowa.v[2]); }

  // force, torque and virial contribution for forces between backbone sites
  delf[0] = 0.0;
  delf[1] = 0.0;
  delf[2] = 0.0;
  delta[0] = 0.0;
  delta[1] = 0.0;
  delta[2] = 0.0;
  deltb[0] = 0.0;
  deltb[1] = 0.0;
  deltb[2] = 0.0;

  // cosphi1 force
  if (cosphi1 != static_cast<KK_FLOAT>(0.0)) {
    finc = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2 * rinv_bkbk);

    delf[0] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[0]*cosphi1 - rowb.v[7]) * finc;
    delf[1] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[1]*cosphi1 - rowb.v[8]) * finc;
    delf[2] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[2]*cosphi1 - rowb.v[9]) * finc;
  }

  // cosphi2 force
  if (cosphi2 != static_cast<KK_FLOAT>(0.0)) {
    finc = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2 * rinv_bkbk);

    delf[0] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[0]*cosphi2 - rowa.v[7]) * finc;
    delf[1] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[1]*cosphi2 - rowa.v[8]) * finc;
    delf[2] += static_cast<KK_ACC_FLOAT>(delr_bkbk_norm[2]*cosphi2 - rowa.v[9]) * finc;
  }

  // accumulate the force and torques of this bond
  fsum[0] += delf[0];
  fsum[1] += delf[1];
  fsum[2] += delf[2];
  tsuma[0] += static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[1];
  tsuma[1] += static_cast<KK_ACC_FLOAT>(ra_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[2];
  tsuma[2] += static_cast<KK_ACC_FLOAT>(ra_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(ra_cbk[1])*delf[0];
  tsumb[0] += static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[2] - static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[1];
  tsumb[1] += static_cast<KK_ACC_FLOAT>(rb_cbk[2])*delf[0] - static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[2];
  tsumb[2] += static_cast<KK_ACC_FLOAT>(rb_cbk[0])*delf[1] - static_cast<KK_ACC_FLOAT>(rb_cbk[1])*delf[0];

  // increment viral only
  if (EVFLAG) { ev_tally_xyz(ev, a, b, nlocal, NEWTON_BOND, 0.0, delf[0], delf[1], delf[2], \
    rowb.v[0]-rowa.v[0], rowb.v[1]-rowa.v[1], rowb.v[2]-rowa.v[2]); }

  // pure torques not expressible as r x f

  delta[0] = 0.0;
  delta[1] = 0.0;
  delta[2] = 0.0;
  deltb[0] = 0.0;
  deltb[1] = 0.0;
  deltb[2] = 0.0;

  // theta4 torque
  if (theta4 != static_cast<KK_FLOAT>(0.0)) {
    tpair = static_cast<KK_ACC_FLOAT>(-f1 * df4t4 * f4t5 * f4t6 * f5c1 * f5c2);
    t4dir[0] = rowa.v[11] * rowb.v[12] - rowa.v[12] * rowb.v[11];
    t4dir[1] = rowa.v[12] * rowb.v[10] - rowa.v[10] * rowb.v[12];
    t4dir[2] = rowa.v[10] * rowb.v[11] - rowa.v[11] * rowb.v[10];
    delta[0] += static_cast<KK_ACC_FLOAT>(t4dir[0]) * tpair;
    delta[1] += static_cast<KK_ACC_FLOAT>(t4dir[1]) * tpair;
    delta[2] += static_cast<KK_ACC_FLOAT>(t4dir[2]) * tpair;
    deltb[0] += static_cast<KK_ACC_FLOAT>(t4dir[0]) * tpair;
    deltb[1] += static_cast<KK_ACC_FLOAT>(t4dir[1]) * tpair;
    deltb[2] += static_cast<KK_ACC_FLOAT>(t4dir[2]) * tpair;
  }

  // theta5p torque
  if (theta5p != static_cast<KK_FLOAT>(0.0)) {
    tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * df4t5 * f4t6 * f5c1 * f5c2);
    t5pdir[0] = delr_stkstk_norm[1] * rowb.v[12] - delr_stkstk_norm[2] * rowb.v[11];
    t5pdir[1] = delr_stkstk_norm[2] * rowb.v[10] - delr_stkstk_norm[0] * rowb.v[12];
    t5pdir[2] = delr_stkstk_norm[0] * rowb.v[11] - delr_stkstk_norm[1] * rowb.v[10];
    deltb[0] += static_cast<KK_ACC_FLOAT>(t5pdir[0]) * tpair;
    deltb[1] += static_cast<KK_ACC_FLOAT>(t5pdir[1]) * tpair;
    deltb[2] += static_cast<KK_ACC_FLOAT>(t5pdir[2]) * tpair;
  }

  // theta6p torque
  if (theta6p != static_cast<KK_FLOAT>(0.0)) {
    tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * df4t6 * f5c1 * f5c2);
    t6pdir[0] = delr_stkstk_norm[1] * rowa.v[12] - delr_stkstk_norm[2] * rowa.v[11];
    t6pdir[1] = delr_stkstk_norm[2] * rowa.v[10] - delr_stkstk_norm[0] * rowa.v[12];
    t6pdir[2] = delr_stkstk_norm[0] * rowa.v[11] - delr_stkstk_norm[1] * rowa.v[10];
    delta[0] -= static_cast<KK_ACC_FLOAT>(t6pdir[0]) * tpair;
    delta[1] -= static_cast<KK_ACC_FLOAT>(t6pdir[1]) * tpair;
    delta[2] -= static_cast<KK_ACC_FLOAT>(t6pdir[2]) * tpair;
  }

  // cosphi1 torque
  if (cosphi1 != static_cast<KK_FLOAT>(0.0)) {
    tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * f4t6 * df5c1 * f5c2);
    cosphi1dir[0] = delr_bkbk_norm[1] * rowb.v[9] - delr_bkbk_norm[2] * rowb.v[8];
    cosphi1dir[1] = delr_bkbk_norm[2] * rowb.v[7] - delr_bkbk_norm[0] * rowb.v[9];
    cosphi1dir[2] = delr_bkbk_norm[0] * rowb.v[8] - delr_bkbk_norm[1] * rowb.v[7];
    deltb[0] += static_cast<KK_ACC_FLOAT>(cosphi1dir[0]) * tpair;
    deltb[1] += static_cast<KK_ACC_FLOAT>(cosphi1dir[1]) * tpair;
    deltb[2] += static_cast<KK_ACC_FLOAT>(cosphi1dir[2]) * tpair;
  }

  // cosphi2 torque
  if (cosphi2 != static_cast<KK_FLOAT>(0.0)) {
    tpair = static_cast<KK_ACC_FLOAT>(-f1 * f4t4 * f4t5 * f4t6 * f5c1 * df5c2);
    cosphi2dir[0] = delr_bkbk_norm[1] * rowa.v[9] - delr_bkbk_norm[2] * rowa.v[8];
    cosphi2dir[1] = delr_bkbk_norm[2] * rowa.v[7] - delr_bkbk_norm[0] * rowa.v[9];
    cosphi2dir[2] = delr_bkbk_norm[0] * rowa.v[8] - delr_bkbk_norm[1] * rowa.v[7];
    delta[0] -= static_cast<KK_ACC_FLOAT>(cosphi2dir[0]) * tpair;
    delta[1] -= static_cast<KK_ACC_FLOAT>(cosphi2dir[1]) * tpair;
    delta[2] -= static_cast<KK_ACC_FLOAT>(cosphi2dir[2]) * tpair;
  }

  // increment forces and torques, once per atom and component
  if ( NEWTON_BOND || a < nlocal ) {
    Kokkos::atomic_add(&f(a,0), -fsum[0]);
    Kokkos::atomic_add(&f(a,1), -fsum[1]);
    Kokkos::atomic_add(&f(a,2), -fsum[2]);
    Kokkos::atomic_add(&torque(a,0), -(tsuma[0] + delta[0]));
    Kokkos::atomic_add(&torque(a,1), -(tsuma[1] + delta[1]));
    Kokkos::atomic_add(&torque(a,2), -(tsuma[2] + delta[2]));
  }
  if ( NEWTON_BOND || b < nlocal ) {
    Kokkos::atomic_add(&f(b,0), fsum[0]);
    Kokkos::atomic_add(&f(b,1), fsum[1]);
    Kokkos::atomic_add(&f(b,2), fsum[2]);
    Kokkos::atomic_add(&torque(b,0), tsumb[0] + deltb[0]);
    Kokkos::atomic_add(&torque(b,1), tsumb[1] + deltb[1]);
    Kokkos::atomic_add(&torque(b,2), tsumb[2] + deltb[2]);
  }
}

template<class DeviceType>
template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION
void PairOxdnaStkKokkos<DeviceType>::operator()(TagPairOxdnaStkCompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, \
  const int &in) const
{
  EV_FLOAT ev;
  this->template operator()<OXDNAFLAG,NEWTON_BOND,EVFLAG>\
  (TagPairOxdnaStkCompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>(),in,ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::allocate()
{
  PairOxdnaStk::allocate();

  int n = atom->ntypes;

  k_params2_st = decltype(k_params2_st)("PairOxdnaStkKokkos:params2_st", n+1, n+1);
  k_params4_st = decltype(k_params4_st)("PairOxdnaStkKokkos:params4_st", n+1, n+1, n+1, n+1);





  d_params2_st = k_params2_st.template view<DeviceType>();
  d_params4_st = k_params4_st.template view<DeviceType>();




}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::settings(int narg, char **/*arg*/)
{
  if (narg != 0) error->all(FLERR, "The oxDNA and oxRNA pair styles do not take any arguments");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::init_style()
{
  // the internal helper fixes are always created for the default KOKKOS variant,
  // so /kk/host styles cannot work with them when LAMMPS is compiled for a GPU

  if (std::is_same_v<DeviceType, LMPHostType> && !std::is_same_v<DeviceType, LMPDeviceType>)
    error->all(FLERR, "The /kk/host variants of the CG-DNA styles are not supported "
               "when LAMMPS is compiled for a GPU");

  // atoms may have been reordered since the last run, so force a rebuild
  // of the cached prime neighbor table in the next compute()

  last_prime_neighs_bond_nbuild = -1;

  // no neighbor list: the stacking interaction loops over the bond list

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0) error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);

  auto prime_fixes = modify->get_fix_by_style("^OXDNA/PRIME_NEIGHS/kk");
  if (prime_fixes.size() == 0)
    fix_oxdna_prime_neighsKK =
      dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(modify->add_fix("prime_neighs_kk all OXDNA/PRIME_NEIGHS/kk"));
  else
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(prime_fixes[0]);

  if (!fix_oxdna_prime_neighsKK)
    error->all(FLERR, "Fix OXDNA/PRIME_NEIGHS/kk not found");

}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxdnaStkKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxdnaStk::init_one(i,j);

  // All non-tetramer Kokkos views are set here within ::init_one, and
  // the tetramer Kokkos views are set within ::coeff

  // Assign directionally: [i][j] gets [i][j], [j][i] gets [j][i]
  k_params2_st.view_host()(i,j).epsilon_st = static_cast<KK_FLOAT>(epsilon_st[i][j]);
  k_params2_st.view_host()(j,i).epsilon_st = static_cast<KK_FLOAT>(epsilon_st[j][i]);
  k_params2_st.view_host()(i,j).a_st = static_cast<KK_FLOAT>(a_st[i][j]);
  k_params2_st.view_host()(j,i).a_st = static_cast<KK_FLOAT>(a_st[j][i]);
  k_params2_st.view_host()(i,j).b_st_lo = static_cast<KK_FLOAT>(b_st_lo[i][j]);
  k_params2_st.view_host()(j,i).b_st_lo = static_cast<KK_FLOAT>(b_st_lo[j][i]);
  k_params2_st.view_host()(i,j).b_st_hi = static_cast<KK_FLOAT>(b_st_hi[i][j]);
  k_params2_st.view_host()(j,i).b_st_hi = static_cast<KK_FLOAT>(b_st_hi[j][i]);

  k_params2_st.view_host()(i,j).theta_st4_0 = static_cast<KK_FLOAT>(theta_st4_0[i][j]);
  k_params2_st.view_host()(j,i).theta_st4_0 = static_cast<KK_FLOAT>(theta_st4_0[j][i]);

  k_params2_st.view_host()(i,j).a_st5 = static_cast<KK_FLOAT>(a_st5[i][j]);
  k_params2_st.view_host()(j,i).a_st5 = static_cast<KK_FLOAT>(a_st5[j][i]);
  k_params2_st.view_host()(i,j).theta_st5_0 = static_cast<KK_FLOAT>(theta_st5_0[i][j]);
  k_params2_st.view_host()(j,i).theta_st5_0 = static_cast<KK_FLOAT>(theta_st5_0[j][i]);
  k_params2_st.view_host()(i,j).dtheta_st5_ast = static_cast<KK_FLOAT>(dtheta_st5_ast[i][j]);
  k_params2_st.view_host()(j,i).dtheta_st5_ast = static_cast<KK_FLOAT>(dtheta_st5_ast[j][i]);
  k_params2_st.view_host()(i,j).b_st5 = static_cast<KK_FLOAT>(b_st5[i][j]);
  k_params2_st.view_host()(j,i).b_st5 = static_cast<KK_FLOAT>(b_st5[j][i]);
  k_params2_st.view_host()(i,j).dtheta_st5_c = static_cast<KK_FLOAT>(dtheta_st5_c[i][j]);
  k_params2_st.view_host()(j,i).dtheta_st5_c = static_cast<KK_FLOAT>(dtheta_st5_c[j][i]);

  k_params2_st.view_host()(i,j).a_st6 = static_cast<KK_FLOAT>(a_st6[i][j]);
  k_params2_st.view_host()(j,i).a_st6 = static_cast<KK_FLOAT>(a_st6[j][i]);
  k_params2_st.view_host()(i,j).theta_st6_0 = static_cast<KK_FLOAT>(theta_st6_0[i][j]);
  k_params2_st.view_host()(j,i).theta_st6_0 = static_cast<KK_FLOAT>(theta_st6_0[j][i]);
  k_params2_st.view_host()(i,j).dtheta_st6_ast = static_cast<KK_FLOAT>(dtheta_st6_ast[i][j]);
  k_params2_st.view_host()(j,i).dtheta_st6_ast = static_cast<KK_FLOAT>(dtheta_st6_ast[j][i]);
  k_params2_st.view_host()(i,j).b_st6 = static_cast<KK_FLOAT>(b_st6[i][j]);
  k_params2_st.view_host()(j,i).b_st6 = static_cast<KK_FLOAT>(b_st6[j][i]);
  k_params2_st.view_host()(i,j).dtheta_st6_c = static_cast<KK_FLOAT>(dtheta_st6_c[i][j]);
  k_params2_st.view_host()(j,i).dtheta_st6_c = static_cast<KK_FLOAT>(dtheta_st6_c[j][i]);

  k_params2_st.view_host()(i,j).a_st1 = static_cast<KK_FLOAT>(a_st1[i][j]);
  k_params2_st.view_host()(j,i).a_st1 = static_cast<KK_FLOAT>(a_st1[j][i]);
  k_params2_st.view_host()(i,j).cosphi_st1_ast = static_cast<KK_FLOAT>(cosphi_st1_ast[i][j]);
  k_params2_st.view_host()(j,i).cosphi_st1_ast = static_cast<KK_FLOAT>(cosphi_st1_ast[j][i]);
  k_params2_st.view_host()(i,j).b_st1 = static_cast<KK_FLOAT>(b_st1[i][j]);
  k_params2_st.view_host()(j,i).b_st1 = static_cast<KK_FLOAT>(b_st1[j][i]);
  k_params2_st.view_host()(i,j).cosphi_st1_c = static_cast<KK_FLOAT>(cosphi_st1_c[i][j]);
  k_params2_st.view_host()(j,i).cosphi_st1_c = static_cast<KK_FLOAT>(cosphi_st1_c[j][i]);
  k_params2_st.view_host()(i,j).a_st2 = static_cast<KK_FLOAT>(a_st2[i][j]);
  k_params2_st.view_host()(j,i).a_st2 = static_cast<KK_FLOAT>(a_st2[j][i]);
  k_params2_st.view_host()(i,j).cosphi_st2_ast = static_cast<KK_FLOAT>(cosphi_st2_ast[i][j]);
  k_params2_st.view_host()(j,i).cosphi_st2_ast = static_cast<KK_FLOAT>(cosphi_st2_ast[j][i]);
  k_params2_st.view_host()(i,j).b_st2 = static_cast<KK_FLOAT>(b_st2[i][j]);
  k_params2_st.view_host()(j,i).b_st2 = static_cast<KK_FLOAT>(b_st2[j][i]);
  k_params2_st.view_host()(i,j).cosphi_st2_c = static_cast<KK_FLOAT>(cosphi_st2_c[i][j]);
  k_params2_st.view_host()(j,i).cosphi_st2_c = static_cast<KK_FLOAT>(cosphi_st2_c[j][i]);

  k_params2_st.modify_host();





  // Sync to device
  k_params2_st.template sync<DeviceType>();





  // "cutone" is max of "cut_st_hc[a][i][j][b]", sets the master list distance cutoff
  return cutone;
}

/* ----------------------------------------------------------------------
   Helper function to set the tetramer Kokkos views within ::coeff
   Is used within child stk/kk classes too.
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::coeff_set_tetramers_kokkos(int narg, char **arg)
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
          k_params4_st.view_host()(i,j,k,l).cut_st_0 = static_cast<KK_FLOAT>(cut_st_0[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cut_st_c = static_cast<KK_FLOAT>(cut_st_c[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cut_st_lo = static_cast<KK_FLOAT>(cut_st_lo[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cut_st_hi = static_cast<KK_FLOAT>(cut_st_hi[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cut_st_lc = static_cast<KK_FLOAT>(cut_st_lc[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cut_st_hc = static_cast<KK_FLOAT>(cut_st_hc[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).shift_st = static_cast<KK_FLOAT>(shift_st[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).cutsq_st_hc = static_cast<KK_FLOAT>(cutsq_st_hc[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).a_st4 = static_cast<KK_FLOAT>(a_st4[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).dtheta_st4_ast = static_cast<KK_FLOAT>(dtheta_st4_ast[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).b_st4 = static_cast<KK_FLOAT>(b_st4[i][j][k][l]);
          k_params4_st.view_host()(i,j,k,l).dtheta_st4_c = static_cast<KK_FLOAT>(dtheta_st4_c[i][j][k][l]);
        }
      }
    }
  }

  k_params4_st.modify_host();

  // Sync to device
  k_params4_st.template sync<DeviceType>();

  // check whether the tetramer-indexed coefficients used by the compute kernel
  // depend on the 3'/5' context types (first and last index) at all

  const int ntypes = atom->ntypes;
  tetramer_uniform = 1;
  for (int j = 1; j <= ntypes && tetramer_uniform; j++) {
    for (int k = 1; k <= ntypes && tetramer_uniform; k++) {
      for (int i = 0; i <= ntypes && tetramer_uniform; i++) {
        for (int l = 0; l <= ntypes; l++) {
          if ((k_params4_st.view_host()(i,j,k,l).cut_st_0 != k_params4_st.view_host()(0,j,k,0).cut_st_0) ||
              (k_params4_st.view_host()(i,j,k,l).cut_st_lc != k_params4_st.view_host()(0,j,k,0).cut_st_lc) ||
              (k_params4_st.view_host()(i,j,k,l).cut_st_hc != k_params4_st.view_host()(0,j,k,0).cut_st_hc) ||
              (k_params4_st.view_host()(i,j,k,l).cut_st_lo != k_params4_st.view_host()(0,j,k,0).cut_st_lo) ||
              (k_params4_st.view_host()(i,j,k,l).cut_st_hi != k_params4_st.view_host()(0,j,k,0).cut_st_hi) ||
              (k_params4_st.view_host()(i,j,k,l).shift_st != k_params4_st.view_host()(0,j,k,0).shift_st) ||
              (k_params4_st.view_host()(i,j,k,l).a_st4 != k_params4_st.view_host()(0,j,k,0).a_st4) ||
              (k_params4_st.view_host()(i,j,k,l).b_st4 != k_params4_st.view_host()(0,j,k,0).b_st4) ||
              (k_params4_st.view_host()(i,j,k,l).dtheta_st4_ast != k_params4_st.view_host()(0,j,k,0).dtheta_st4_ast) ||
              (k_params4_st.view_host()(i,j,k,l).dtheta_st4_c != k_params4_st.view_host()(0,j,k,0).dtheta_st4_c)) {
            tetramer_uniform = 0;
            break;
          }
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   The tetramer Kokkos views are set here within ::coeff, and the
   non-tetramer Kokkos views are set within ::init_one
------------------------------------------------------------------------- */

template<class DeviceType>
void PairOxdnaStkKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairOxdnaStk::coeff(narg,arg);

  coeff_set_tetramers_kokkos(narg,arg);
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators

   NOTE: Although this is a pair style interaction, the algorithm below
   follows the virial incrementation of the bond style. This is because
   the bond topology is used in the main compute loop.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void PairOxdnaStkKokkos<DeviceType>::ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,\
      const int &nlocal, const int &newton_bond, const KK_FLOAT &evdwl,\
      const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,\
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const
{
  KK_ACC_FLOAT evdwlhalf;
  KK_ACC_FLOAT v[6];

  // The eatom and vatom arrays are atomic
  Kokkos::View<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout, \
      typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout, \
      typename KKDevice<DeviceType>::value,Kokkos::MemoryTraits<Kokkos::Atomic|Kokkos::Unmanaged> > v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond) ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      else {
        evdwlhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*evdwl);
        if (i < nlocal) ev.evdwl += evdwlhalf;
        if (j < nlocal) ev.evdwl += evdwlhalf;
      }
    }
    if (eflag_atom) {
      evdwlhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5)*evdwl);
      if (newton_bond || i < nlocal) v_eatom[i] += evdwlhalf;
      if (newton_bond || j < nlocal) v_eatom[j] += evdwlhalf;
    }
  }

  if (vflag_either) {
    v[0] = static_cast<KK_ACC_FLOAT>(delx) * fx;
    v[1] = static_cast<KK_ACC_FLOAT>(dely) * fy;
    v[2] = static_cast<KK_ACC_FLOAT>(delz) * fz;
    v[3] = static_cast<KK_ACC_FLOAT>(delx) * fy;
    v[4] = static_cast<KK_ACC_FLOAT>(delx) * fz;
    v[5] = static_cast<KK_ACC_FLOAT>(dely) * fz;

    if (vflag_global) {
      if (newton_bond) {
        ev.v[0] += v[0];
        ev.v[1] += v[1];
        ev.v[2] += v[2];
        ev.v[3] += v[3];
        ev.v[4] += v[4];
        ev.v[5] += v[5];
      } else {
        if (i < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v[0];
          ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v[1];
          ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v[2];
          ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v[3];
          ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v[4];
          ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v[5];
        }
        if (j < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5)*v[0];
          ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5)*v[1];
          ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5)*v[2];
          ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5)*v[3];
          ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5)*v[4];
          ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5)*v[5];
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i,0) += static_cast<KK_ACC_FLOAT>(0.5)*v[0];
        v_vatom(i,1) += static_cast<KK_ACC_FLOAT>(0.5)*v[1];
        v_vatom(i,2) += static_cast<KK_ACC_FLOAT>(0.5)*v[2];
        v_vatom(i,3) += static_cast<KK_ACC_FLOAT>(0.5)*v[3];
        v_vatom(i,4) += static_cast<KK_ACC_FLOAT>(0.5)*v[4];
        v_vatom(i,5) += static_cast<KK_ACC_FLOAT>(0.5)*v[5];
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j,0) += static_cast<KK_ACC_FLOAT>(0.5)*v[0];
        v_vatom(j,1) += static_cast<KK_ACC_FLOAT>(0.5)*v[1];
        v_vatom(j,2) += static_cast<KK_ACC_FLOAT>(0.5)*v[2];
        v_vatom(j,3) += static_cast<KK_ACC_FLOAT>(0.5)*v[3];
        v_vatom(j,4) += static_cast<KK_ACC_FLOAT>(0.5)*v[4];
        v_vatom(j,5) += static_cast<KK_ACC_FLOAT>(0.5)*v[5];
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairOxdnaStkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxdnaStkKokkos<LMPHostType>;
#endif
}
