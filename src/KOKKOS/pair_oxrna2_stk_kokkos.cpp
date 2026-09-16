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

#include "pair_oxrna2_stk_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "memory_kokkos.h"
#include "modify.h"
#include "neigh_request.h"
#include "neighbor_kokkos.h"

#include "fix_oxdna_lrf_kokkos.h"
#include "fix_oxdna_prime_neighs_kokkos.h"
#include "mf_oxdna_kokkos.h"

using namespace LAMMPS_NS;
using namespace MFOxdnaKokkos;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxrna2StkKokkos<DeviceType>::PairOxrna2StkKokkos(LAMMPS *lmp) : PairOxrna2Stk(lmp)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  neighborKK = (NeighborKokkos *) neighbor;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  // Internal FixOxdnaLRFKokkos already syncs all read masks that do not
  // change between pair/bond styles.
  datamask_read = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | TORQUE_MASK | ENERGY_MASK | VIRIAL_MASK;

  fix_oxdna_lrfKK = nullptr;
  fix_oxdna_prime_neighsKK = nullptr;
  last_prime_neighs_bond_lastcall = -1;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairOxrna2StkKokkos<DeviceType>::~PairOxrna2StkKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom, eatom);
    memoryKK->destroy_kokkos(k_vatom, vatom);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2StkKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  ev_init(eflag, vflag, 0);

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

  atomKK->sync(execution_space, datamask_read);

  if (eflag || vflag)
    atomKK->modified(execution_space, datamask_modify);
  else
    atomKK->modified(execution_space, F_MASK | TORQUE_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  torque = atomKK->k_torque.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();

  nlocal = atom->nlocal;
  newton_bond = force->newton_bond;
  nbondlist = neighborKK->nbondlist;

  // Keep bond-context precompute aligned with the current neighbor-list epoch.
  if (last_prime_neighs_bond_lastcall != neighbor->lastcall) {
    fix_oxdna_prime_neighsKK->compute_prime_neighs_bond();
    last_prime_neighs_bond_lastcall = neighbor->lastcall;
  }

  d_prime_neighs_bond = fix_oxdna_prime_neighsKK->d_prime_neighs_bond;

  copymode = 1;

  // d_n(x/y/z)_xtrct = extracted local unit vectors in lab frame from fix_oxdna_lrf_kokkos.
  d_nx_xtrct = fix_oxdna_lrfKK->k_nx.template view<DeviceType>();
  d_ny_xtrct = fix_oxdna_lrfKK->k_ny.template view<DeviceType>();
  d_nz_xtrct = fix_oxdna_lrfKK->k_nz.template view<DeviceType>();

  EV_FLOAT ev;

  if (evflag) {
    if (newton_bond) {
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType, TagPairOxrna2StkCompute<1, 1>>(0, nbondlist), *this,
          ev);
    } else {
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<DeviceType, TagPairOxrna2StkCompute<0, 1>>(0, nbondlist), *this,
          ev);
    }
  } else {
    if (newton_bond) {
      Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceType, TagPairOxrna2StkCompute<1, 0>>(0, nbondlist), *this);
    } else {
      Kokkos::parallel_for(
          Kokkos::RangePolicy<DeviceType, TagPairOxrna2StkCompute<0, 0>>(0, nbondlist), *this);
    }
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
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  copymode = 0;
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairOxrna2StkKokkos<DeviceType>::operator()(
    TagPairOxrna2StkCompute<NEWTON_BOND, EVFLAG>, const int &in, EV_FLOAT &ev) const
{
  // The f and torque arrays are atomic.
  Kokkos::View<KK_ACC_FLOAT *[3], typename DAT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      a_f = f;
  Kokkos::View<KK_ACC_FLOAT *[3], typename DAT::t_kkacc_1d_3::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      a_torque = torque;

  // Precomputed ordering guarantees a=3' and b=5'.
  const int a = d_prime_neighs_bond(in, 0);
  const int b = d_prime_neighs_bond(in, 1);
  const int atype = type(a);
  const int btype = type(b);

  constexpr KK_FLOAT dx_cbk_oxrna2 = -0.4;
  constexpr KK_FLOAT dz_cbk_oxrna2 = +0.2;
  constexpr KK_FLOAT dx_cstk_3p_oxrna2 = +0.4;
  constexpr KK_FLOAT dy_cstk_3p_oxrna2 = +0.1;
  constexpr KK_FLOAT dx_cstk_5p_oxrna2 = +0.124906078525;
  constexpr KK_FLOAT dy_cstk_5p_oxrna2 = -0.00866274917473;

  constexpr KK_FLOAT d3p_x = -0.462510;
  constexpr KK_FLOAT d3p_y = -0.528218;
  constexpr KK_FLOAT d3p_z = +0.712089;
  constexpr KK_FLOAT d5p_x = -0.104402;
  constexpr KK_FLOAT d5p_y = -0.841783;
  constexpr KK_FLOAT d5p_z = +0.529624;

  KK_FLOAT ra_cstk[3], rb_cstk[3];
  KK_FLOAT ra_cbk[3], rb_cbk[3];
  KK_FLOAT ax[3], ay[3], az[3];
  KK_FLOAT bx[3], by[3], bz[3];
  KK_FLOAT aux3p[3], aux5p[3];

  KK_ACC_FLOAT delf[3], delta[3], deltb[3];
  KK_ACC_FLOAT evdwl, finc, tpair;

  KK_FLOAT delr_bkbk[3], delr_bkbk_norm[3], rsq_bkbk, r_bkbk, rinv_bkbk;
  KK_FLOAT delr_stkstk[3], delr_stkstk_norm[3], rsq_stkstk, r_stkstk, rinv_stkstk;
  KK_FLOAT theta5p, t5pdir[3], cost5p;
  KK_FLOAT theta6p, t6pdir[3], cost6p;
  KK_FLOAT theta9, t9dir[3], cost9;
  KK_FLOAT theta10, t10dir[3], cost10;
  KK_FLOAT cosphi1, cosphi2, cosphi1dir[3], cosphi2dir[3];

  KK_FLOAT f1, f4t5, f4t6, f4t9, f4t10, f5c1, f5c2;
  KK_FLOAT df1, df4t5, df4t6, df4t9, df4t10, df5c1, df5c2;

  ax[0] = d_nx_xtrct(a, 0);
  ax[1] = d_nx_xtrct(a, 1);
  ax[2] = d_nx_xtrct(a, 2);
  ay[0] = d_ny_xtrct(a, 0);
  ay[1] = d_ny_xtrct(a, 1);
  ay[2] = d_ny_xtrct(a, 2);
  az[0] = d_nz_xtrct(a, 0);
  az[1] = d_nz_xtrct(a, 1);
  az[2] = d_nz_xtrct(a, 2);

  bx[0] = d_nx_xtrct(b, 0);
  bx[1] = d_nx_xtrct(b, 1);
  bx[2] = d_nx_xtrct(b, 2);
  by[0] = d_ny_xtrct(b, 0);
  by[1] = d_ny_xtrct(b, 1);
  by[2] = d_ny_xtrct(b, 2);
  bz[0] = d_nz_xtrct(b, 0);
  bz[1] = d_nz_xtrct(b, 1);
  bz[2] = d_nz_xtrct(b, 2);

  // vector COM a - 5'-stacking site a
  ra_cstk[0] = dx_cstk_5p_oxrna2 * ax[0] + dy_cstk_5p_oxrna2 * ay[0];
  ra_cstk[1] = dx_cstk_5p_oxrna2 * ax[1] + dy_cstk_5p_oxrna2 * ay[1];
  ra_cstk[2] = dx_cstk_5p_oxrna2 * ax[2] + dy_cstk_5p_oxrna2 * ay[2];

  // vector COM b - 3'-stacking site b
  rb_cstk[0] = dx_cstk_3p_oxrna2 * bx[0] + dy_cstk_3p_oxrna2 * by[0];
  rb_cstk[1] = dx_cstk_3p_oxrna2 * bx[1] + dy_cstk_3p_oxrna2 * by[1];
  rb_cstk[2] = dx_cstk_3p_oxrna2 * bx[2] + dy_cstk_3p_oxrna2 * by[2];

  // vector 5'-stacking site a to 3'-stacking site b
  delr_stkstk[0] = x(b, 0) + rb_cstk[0] - x(a, 0) - ra_cstk[0];
  delr_stkstk[1] = x(b, 1) + rb_cstk[1] - x(a, 1) - ra_cstk[1];
  delr_stkstk[2] = x(b, 2) + rb_cstk[2] - x(a, 2) - ra_cstk[2];

  rsq_stkstk = Kokkos::fma(delr_stkstk[0], delr_stkstk[0], Kokkos::fma(delr_stkstk[1], delr_stkstk[1], delr_stkstk[2]*delr_stkstk[2]));
  r_stkstk = sqrtf(rsq_stkstk);
  rinv_stkstk = 1.0 / r_stkstk;

  delr_stkstk_norm[0] = delr_stkstk[0] * rinv_stkstk;
  delr_stkstk_norm[1] = delr_stkstk[1] * rinv_stkstk;
  delr_stkstk_norm[2] = delr_stkstk[2] * rinv_stkstk;

  // vector COM a - backbone site a
  ra_cbk[0] = dx_cbk_oxrna2 * ax[0] + dz_cbk_oxrna2 * az[0];
  ra_cbk[1] = dx_cbk_oxrna2 * ax[1] + dz_cbk_oxrna2 * az[1];
  ra_cbk[2] = dx_cbk_oxrna2 * ax[2] + dz_cbk_oxrna2 * az[2];

  // vector COM b - backbone site b
  rb_cbk[0] = dx_cbk_oxrna2 * bx[0] + dz_cbk_oxrna2 * bz[0];
  rb_cbk[1] = dx_cbk_oxrna2 * bx[1] + dz_cbk_oxrna2 * bz[1];
  rb_cbk[2] = dx_cbk_oxrna2 * bx[2] + dz_cbk_oxrna2 * bz[2];

  // vector backbone site a to b
  delr_bkbk[0] = x(b, 0) + rb_cbk[0] - x(a, 0) - ra_cbk[0];
  delr_bkbk[1] = x(b, 1) + rb_cbk[1] - x(a, 1) - ra_cbk[1];
  delr_bkbk[2] = x(b, 2) + rb_cbk[2] - x(a, 2) - ra_cbk[2];

  rsq_bkbk = Kokkos::fma(delr_bkbk[0], delr_bkbk[0], Kokkos::fma(delr_bkbk[1], delr_bkbk[1], delr_bkbk[2]*delr_bkbk[2]));
  r_bkbk = sqrtf(rsq_bkbk);
  rinv_bkbk = 1.0 / r_bkbk;

  delr_bkbk_norm[0] = delr_bkbk[0] * rinv_bkbk;
  delr_bkbk_norm[1] = delr_bkbk[1] * rinv_bkbk;
  delr_bkbk_norm[2] = delr_bkbk[2] * rinv_bkbk;

  f1 = F1_KK(r_stkstk, d_epsilon_st(atype, btype), d_a_st(atype, btype), d_cut_st_0(atype, btype),
             d_cut_st_lc(atype, btype), d_cut_st_hc(atype, btype), d_cut_st_lo(atype, btype),
             d_cut_st_hi(atype, btype), d_b_st_lo(atype, btype), d_b_st_hi(atype, btype),
             d_shift_st(atype, btype));

  if (f1 != 0.0) {
    cost5p = delr_stkstk_norm[0] * bz[0] + delr_stkstk_norm[1] * bz[1] +
        delr_stkstk_norm[2] * bz[2];
    if (cost5p > 1.0) cost5p = 1.0;
    if (cost5p < -1.0) cost5p = -1.0;
    theta5p = acos(cost5p);

    f4t5 = F4_KK(theta5p, d_a_st5(atype, btype), d_theta_st5_0(atype, btype),
                 d_dtheta_st5_ast(atype, btype), d_b_st5(atype, btype),
                 d_dtheta_st5_c(atype, btype));

    if (f4t5 != 0.0) {
      cost6p = delr_stkstk_norm[0] * az[0] + delr_stkstk_norm[1] * az[1] +
          delr_stkstk_norm[2] * az[2];
      if (cost6p > 1.0) cost6p = 1.0;
      if (cost6p < -1.0) cost6p = -1.0;
      theta6p = acos(cost6p);

      aux5p[0] = d5p_x * ax[0] + d5p_y * ay[0] + d5p_z * az[0];
      aux5p[1] = d5p_x * ax[1] + d5p_y * ay[1] + d5p_z * az[1];
      aux5p[2] = d5p_x * ax[2] + d5p_y * ay[2] + d5p_z * az[2];

      aux3p[0] = d3p_x * bx[0] + d3p_y * by[0] + d3p_z * bz[0];
      aux3p[1] = d3p_x * bx[1] + d3p_y * by[1] + d3p_z * bz[1];
      aux3p[2] = d3p_x * bx[2] + d3p_y * by[2] + d3p_z * bz[2];

      cost9 = delr_bkbk_norm[0] * aux3p[0] + delr_bkbk_norm[1] * aux3p[1] +
          delr_bkbk_norm[2] * aux3p[2];
      if (cost9 > 1.0) cost9 = 1.0;
      if (cost9 < -1.0) cost9 = -1.0;
      theta9 = acos(cost9);

      cost10 = delr_bkbk_norm[0] * aux5p[0] + delr_bkbk_norm[1] * aux5p[1] +
          delr_bkbk_norm[2] * aux5p[2];
      if (cost10 > 1.0) cost10 = 1.0;
      if (cost10 < -1.0) cost10 = -1.0;
      theta10 = acos(cost10);

      cosphi1 = delr_bkbk_norm[0] * by[0] + delr_bkbk_norm[1] * by[1] + delr_bkbk_norm[2] * by[2];
      if (cosphi1 > 1.0) cosphi1 = 1.0;
      if (cosphi1 < -1.0) cosphi1 = -1.0;

      cosphi2 = delr_bkbk_norm[0] * ay[0] + delr_bkbk_norm[1] * ay[1] + delr_bkbk_norm[2] * ay[2];
      if (cosphi2 > 1.0) cosphi2 = 1.0;
      if (cosphi2 < -1.0) cosphi2 = -1.0;

      f4t6 = F4_KK(theta6p, d_a_st6(atype, btype), d_theta_st6_0(atype, btype),
                   d_dtheta_st6_ast(atype, btype), d_b_st6(atype, btype),
                   d_dtheta_st6_c(atype, btype));
      f4t9 = F4_KK(theta9, d_a_st9(atype, btype), d_theta_st9_0(atype, btype),
                   d_dtheta_st9_ast(atype, btype), d_b_st9(atype, btype),
                   d_dtheta_st9_c(atype, btype));
      f4t10 = F4_KK(theta10, d_a_st10(atype, btype), d_theta_st10_0(atype, btype),
                    d_dtheta_st10_ast(atype, btype), d_b_st10(atype, btype),
                    d_dtheta_st10_c(atype, btype));
      f5c1 = F5_KK(-cosphi1, d_a_st1(atype, btype), -d_cosphi_st1_ast(atype, btype),
                   d_b_st1(atype, btype), -d_cosphi_st1_c(atype, btype));
      f5c2 = F5_KK(-cosphi2, d_a_st2(atype, btype), -d_cosphi_st2_ast(atype, btype),
                   d_b_st2(atype, btype), -d_cosphi_st2_c(atype, btype));

      evdwl = f1 * f4t5 * f4t6 * f4t9 * f4t10 * f5c1 * f5c2;

      if (evdwl != 0.0) {
        df1 = DF1_KK(r_stkstk, d_epsilon_st(atype, btype), d_a_st(atype, btype),
                     d_cut_st_0(atype, btype), d_cut_st_lc(atype, btype),
                     d_cut_st_hc(atype, btype), d_cut_st_lo(atype, btype),
                     d_cut_st_hi(atype, btype), d_b_st_lo(atype, btype), d_b_st_hi(atype, btype));
        df4t5 = DF4_KK(theta5p, d_a_st5(atype, btype), d_theta_st5_0(atype, btype),
                       d_dtheta_st5_ast(atype, btype), d_b_st5(atype, btype),
                       d_dtheta_st5_c(atype, btype)) /
            sin(theta5p);
        df4t6 = DF4_KK(theta6p, d_a_st6(atype, btype), d_theta_st6_0(atype, btype),
                       d_dtheta_st6_ast(atype, btype), d_b_st6(atype, btype),
                       d_dtheta_st6_c(atype, btype)) /
            sin(theta6p);
        df4t9 = DF4_KK(theta9, d_a_st9(atype, btype), d_theta_st9_0(atype, btype),
                       d_dtheta_st9_ast(atype, btype), d_b_st9(atype, btype),
                       d_dtheta_st9_c(atype, btype)) /
            sin(theta9);
        df4t10 = DF4_KK(theta10, d_a_st10(atype, btype), d_theta_st10_0(atype, btype),
                        d_dtheta_st10_ast(atype, btype), d_b_st10(atype, btype),
                        d_dtheta_st10_c(atype, btype)) /
            sin(theta10);
        df5c1 = DF5_KK(-cosphi1, d_a_st1(atype, btype), -d_cosphi_st1_ast(atype, btype),
                       d_b_st1(atype, btype), -d_cosphi_st1_c(atype, btype));
        df5c2 = DF5_KK(-cosphi2, d_a_st2(atype, btype), -d_cosphi_st2_ast(atype, btype),
                       d_b_st2(atype, btype), -d_cosphi_st2_c(atype, btype));

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
        finc = -df1 * f4t5 * f4t6 * f4t9 * f4t10 * f5c1 * f5c2;
        delf[0] += delr_stkstk[0] * finc;
        delf[1] += delr_stkstk[1] * finc;
        delf[2] += delr_stkstk[2] * finc;

        // theta5p force
        if (theta5p != 0.0) {
          finc = -f1 * df4t5 * f4t6 * f4t9 * f4t10 * f5c1 * f5c2 * rinv_stkstk;
          delf[0] += (delr_stkstk_norm[0] * cost5p - bz[0]) * finc;
          delf[1] += (delr_stkstk_norm[1] * cost5p - bz[1]) * finc;
          delf[2] += (delr_stkstk_norm[2] * cost5p - bz[2]) * finc;
        }

        // theta6p force
        if (theta6p != 0.0) {
          finc = -f1 * f4t5 * df4t6 * f4t9 * f4t10 * f5c1 * f5c2 * rinv_stkstk;
          delf[0] += (delr_stkstk_norm[0] * cost6p - az[0]) * finc;
          delf[1] += (delr_stkstk_norm[1] * cost6p - az[1]) * finc;
          delf[2] += (delr_stkstk_norm[2] * cost6p - az[2]) * finc;
        }

        if (NEWTON_BOND || a < nlocal) {
          a_f(a, 0) -= delf[0];
          a_f(a, 1) -= delf[1];
          a_f(a, 2) -= delf[2];
          delta[0] = ra_cstk[1] * delf[2] - ra_cstk[2] * delf[1];
          delta[1] = ra_cstk[2] * delf[0] - ra_cstk[0] * delf[2];
          delta[2] = ra_cstk[0] * delf[1] - ra_cstk[1] * delf[0];
          a_torque(a, 0) -= delta[0];
          a_torque(a, 1) -= delta[1];
          a_torque(a, 2) -= delta[2];
        }
        if (NEWTON_BOND || b < nlocal) {
          a_f(b, 0) += delf[0];
          a_f(b, 1) += delf[1];
          a_f(b, 2) += delf[2];
          deltb[0] = rb_cstk[1] * delf[2] - rb_cstk[2] * delf[1];
          deltb[1] = rb_cstk[2] * delf[0] - rb_cstk[0] * delf[2];
          deltb[2] = rb_cstk[0] * delf[1] - rb_cstk[1] * delf[0];
          a_torque(b, 0) += deltb[0];
          a_torque(b, 1) += deltb[1];
          a_torque(b, 2) += deltb[2];
        }

        if (EVFLAG) {
          ev_tally_xyz(ev, b, a, nlocal, NEWTON_BOND, evdwl, delf[0], delf[1], delf[2], x(b, 0) - x(a, 0),
                       x(b, 1) - x(a, 1), x(b, 2) - x(a, 2));
        }

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

        // theta9 force
        if (theta9 != 0.0) {
          finc = -f1 * f4t5 * f4t6 * df4t9 * f4t10 * f5c1 * f5c2 * rinv_bkbk;
          delf[0] += (delr_bkbk_norm[0] * cost9 - aux3p[0]) * finc;
          delf[1] += (delr_bkbk_norm[1] * cost9 - aux3p[1]) * finc;
          delf[2] += (delr_bkbk_norm[2] * cost9 - aux3p[2]) * finc;
        }

        // theta10 force
        if (theta10 != 0.0) {
          finc = -f1 * f4t5 * f4t6 * f4t9 * df4t10 * f5c1 * f5c2 * rinv_bkbk;
          delf[0] += (delr_bkbk_norm[0] * cost10 - aux5p[0]) * finc;
          delf[1] += (delr_bkbk_norm[1] * cost10 - aux5p[1]) * finc;
          delf[2] += (delr_bkbk_norm[2] * cost10 - aux5p[2]) * finc;
        }

        // cosphi1 force
        if (cosphi1 != 0.0) {
          finc = -f1 * f4t5 * f4t6 * f4t9 * f4t10 * df5c1 * f5c2 * rinv_bkbk;
          delf[0] += (delr_bkbk_norm[0] * cosphi1 - by[0]) * finc;
          delf[1] += (delr_bkbk_norm[1] * cosphi1 - by[1]) * finc;
          delf[2] += (delr_bkbk_norm[2] * cosphi1 - by[2]) * finc;
        }

        // cosphi2 force
        if (cosphi2 != 0.0) {
          finc = -f1 * f4t5 * f4t6 * f4t9 * f4t10 * f5c1 * df5c2 * rinv_bkbk;
          delf[0] += (delr_bkbk_norm[0] * cosphi2 - ay[0]) * finc;
          delf[1] += (delr_bkbk_norm[1] * cosphi2 - ay[1]) * finc;
          delf[2] += (delr_bkbk_norm[2] * cosphi2 - ay[2]) * finc;
        }

        if (NEWTON_BOND || a < nlocal) {
          a_f(a, 0) -= delf[0];
          a_f(a, 1) -= delf[1];
          a_f(a, 2) -= delf[2];
          delta[0] = ra_cbk[1] * delf[2] - ra_cbk[2] * delf[1];
          delta[1] = ra_cbk[2] * delf[0] - ra_cbk[0] * delf[2];
          delta[2] = ra_cbk[0] * delf[1] - ra_cbk[1] * delf[0];
          a_torque(a, 0) -= delta[0];
          a_torque(a, 1) -= delta[1];
          a_torque(a, 2) -= delta[2];
        }
        if (NEWTON_BOND || b < nlocal) {
          a_f(b, 0) += delf[0];
          a_f(b, 1) += delf[1];
          a_f(b, 2) += delf[2];
          deltb[0] = rb_cbk[1] * delf[2] - rb_cbk[2] * delf[1];
          deltb[1] = rb_cbk[2] * delf[0] - rb_cbk[0] * delf[2];
          deltb[2] = rb_cbk[0] * delf[1] - rb_cbk[1] * delf[0];
          a_torque(b, 0) += deltb[0];
          a_torque(b, 1) += deltb[1];
          a_torque(b, 2) += deltb[2];
        }

        if (EVFLAG) {
          ev_tally_xyz(ev, b, a, nlocal, NEWTON_BOND, 0.0, delf[0], delf[1], delf[2], x(b, 0) - x(a, 0),
                       x(b, 1) - x(a, 1), x(b, 2) - x(a, 2));
        }

        // pure torques not expressible as r x f
        delta[0] = 0.0;
        delta[1] = 0.0;
        delta[2] = 0.0;
        deltb[0] = 0.0;
        deltb[1] = 0.0;
        deltb[2] = 0.0;

        // theta5p torque
        if (theta5p != 0.0) {
          tpair = -f1 * df4t5 * f4t6 * f4t9 * f4t10 * f5c1 * f5c2;
          t5pdir[0] = delr_stkstk_norm[1] * bz[2] - delr_stkstk_norm[2] * bz[1];
          t5pdir[1] = delr_stkstk_norm[2] * bz[0] - delr_stkstk_norm[0] * bz[2];
          t5pdir[2] = delr_stkstk_norm[0] * bz[1] - delr_stkstk_norm[1] * bz[0];
          deltb[0] += t5pdir[0] * tpair;
          deltb[1] += t5pdir[1] * tpair;
          deltb[2] += t5pdir[2] * tpair;
        }

        // theta6p torque
        if (theta6p != 0.0) {
          tpair = -f1 * f4t5 * df4t6 * f4t9 * f4t10 * f5c1 * f5c2;
          t6pdir[0] = delr_stkstk_norm[1] * az[2] - delr_stkstk_norm[2] * az[1];
          t6pdir[1] = delr_stkstk_norm[2] * az[0] - delr_stkstk_norm[0] * az[2];
          t6pdir[2] = delr_stkstk_norm[0] * az[1] - delr_stkstk_norm[1] * az[0];
          delta[0] -= t6pdir[0] * tpair;
          delta[1] -= t6pdir[1] * tpair;
          delta[2] -= t6pdir[2] * tpair;
        }

        // theta9 torque
        if (theta9 != 0.0) {
          tpair = -f1 * f4t5 * f4t6 * df4t9 * f4t10 * f5c1 * f5c2;
          t9dir[0] = delr_bkbk_norm[1] * aux3p[2] - delr_bkbk_norm[2] * aux3p[1];
          t9dir[1] = delr_bkbk_norm[2] * aux3p[0] - delr_bkbk_norm[0] * aux3p[2];
          t9dir[2] = delr_bkbk_norm[0] * aux3p[1] - delr_bkbk_norm[1] * aux3p[0];
          deltb[0] += t9dir[0] * tpair;
          deltb[1] += t9dir[1] * tpair;
          deltb[2] += t9dir[2] * tpair;
        }

        // theta10 torque
        if (theta10 != 0.0) {
          tpair = -f1 * f4t5 * f4t6 * f4t9 * df4t10 * f5c1 * f5c2;
          t10dir[0] = delr_bkbk_norm[1] * aux5p[2] - delr_bkbk_norm[2] * aux5p[1];
          t10dir[1] = delr_bkbk_norm[2] * aux5p[0] - delr_bkbk_norm[0] * aux5p[2];
          t10dir[2] = delr_bkbk_norm[0] * aux5p[1] - delr_bkbk_norm[1] * aux5p[0];
          delta[0] -= t10dir[0] * tpair;
          delta[1] -= t10dir[1] * tpair;
          delta[2] -= t10dir[2] * tpair;
        }

        // cosphi1 torque
        if (cosphi1 != 0.0) {
          tpair = -f1 * f4t5 * f4t6 * f4t9 * f4t10 * df5c1 * f5c2;
          cosphi1dir[0] = delr_bkbk_norm[1] * by[2] - delr_bkbk_norm[2] * by[1];
          cosphi1dir[1] = delr_bkbk_norm[2] * by[0] - delr_bkbk_norm[0] * by[2];
          cosphi1dir[2] = delr_bkbk_norm[0] * by[1] - delr_bkbk_norm[1] * by[0];
          deltb[0] += cosphi1dir[0] * tpair;
          deltb[1] += cosphi1dir[1] * tpair;
          deltb[2] += cosphi1dir[2] * tpair;
        }

        // cosphi2 torque
        if (cosphi2 != 0.0) {
          tpair = -f1 * f4t5 * f4t6 * f4t9 * f4t10 * f5c1 * df5c2;
          cosphi2dir[0] = delr_bkbk_norm[1] * ay[2] - delr_bkbk_norm[2] * ay[1];
          cosphi2dir[1] = delr_bkbk_norm[2] * ay[0] - delr_bkbk_norm[0] * ay[2];
          cosphi2dir[2] = delr_bkbk_norm[0] * ay[1] - delr_bkbk_norm[1] * ay[0];
          delta[0] -= cosphi2dir[0] * tpair;
          delta[1] -= cosphi2dir[1] * tpair;
          delta[2] -= cosphi2dir[2] * tpair;
        }

        if (NEWTON_BOND || a < nlocal) {
          a_torque(a, 0) -= delta[0];
          a_torque(a, 1) -= delta[1];
          a_torque(a, 2) -= delta[2];
        }
        if (NEWTON_BOND || b < nlocal) {
          a_torque(b, 0) += deltb[0];
          a_torque(b, 1) += deltb[1];
          a_torque(b, 2) += deltb[2];
        }
      }
    }
  }
}

template<class DeviceType>
template<int NEWTON_BOND, int EVFLAG>
KOKKOS_INLINE_FUNCTION void PairOxrna2StkKokkos<DeviceType>::operator()(
    TagPairOxrna2StkCompute<NEWTON_BOND, EVFLAG>, const int &in) const
{
  EV_FLOAT ev;
  this->template operator()<NEWTON_BOND, EVFLAG>(
      TagPairOxrna2StkCompute<NEWTON_BOND, EVFLAG>(), in, ev);
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2StkKokkos<DeviceType>::allocate()
{
  PairOxrna2Stk::allocate();

  int n = atom->ntypes;

  memoryKK->create_kokkos(k_epsilon_st, n + 1, n + 1, "PairOxrna2Stk:epsilon_st");
  memoryKK->create_kokkos(k_a_st, n + 1, n + 1, "PairOxrna2Stk:a_st");
  memoryKK->create_kokkos(k_cut_st_0, n + 1, n + 1, "PairOxrna2Stk:cut_st_0");
  memoryKK->create_kokkos(k_cut_st_c, n + 1, n + 1, "PairOxrna2Stk:cut_st_c");
  memoryKK->create_kokkos(k_cut_st_lo, n + 1, n + 1, "PairOxrna2Stk:cut_st_lo");
  memoryKK->create_kokkos(k_cut_st_hi, n + 1, n + 1, "PairOxrna2Stk:cut_st_hi");
  memoryKK->create_kokkos(k_cut_st_lc, n + 1, n + 1, "PairOxrna2Stk:cut_st_lc");
  memoryKK->create_kokkos(k_cut_st_hc, n + 1, n + 1, "PairOxrna2Stk:cut_st_hc");
  memoryKK->create_kokkos(k_b_st_lo, n + 1, n + 1, "PairOxrna2Stk:b_st_lo");
  memoryKK->create_kokkos(k_b_st_hi, n + 1, n + 1, "PairOxrna2Stk:b_st_hi");
  memoryKK->create_kokkos(k_shift_st, n + 1, n + 1, "PairOxrna2Stk:shift_st");
  memoryKK->create_kokkos(k_cutsq_st_hc, n + 1, n + 1, "PairOxrna2Stk:cutsq_st_hc");

  memoryKK->create_kokkos(k_a_st5, n + 1, n + 1, "PairOxrna2Stk:a_st5");
  memoryKK->create_kokkos(k_theta_st5_0, n + 1, n + 1, "PairOxrna2Stk:theta_st5_0");
  memoryKK->create_kokkos(k_dtheta_st5_ast, n + 1, n + 1, "PairOxrna2Stk:dtheta_st5_ast");
  memoryKK->create_kokkos(k_b_st5, n + 1, n + 1, "PairOxrna2Stk:b_st5");
  memoryKK->create_kokkos(k_dtheta_st5_c, n + 1, n + 1, "PairOxrna2Stk:dtheta_st5_c");

  memoryKK->create_kokkos(k_a_st6, n + 1, n + 1, "PairOxrna2Stk:a_st6");
  memoryKK->create_kokkos(k_theta_st6_0, n + 1, n + 1, "PairOxrna2Stk:theta_st6_0");
  memoryKK->create_kokkos(k_dtheta_st6_ast, n + 1, n + 1, "PairOxrna2Stk:dtheta_st6_ast");
  memoryKK->create_kokkos(k_b_st6, n + 1, n + 1, "PairOxrna2Stk:b_st6");
  memoryKK->create_kokkos(k_dtheta_st6_c, n + 1, n + 1, "PairOxrna2Stk:dtheta_st6_c");

  memoryKK->create_kokkos(k_a_st9, n + 1, n + 1, "PairOxrna2Stk:a_st9");
  memoryKK->create_kokkos(k_theta_st9_0, n + 1, n + 1, "PairOxrna2Stk:theta_st9_0");
  memoryKK->create_kokkos(k_dtheta_st9_ast, n + 1, n + 1, "PairOxrna2Stk:dtheta_st9_ast");
  memoryKK->create_kokkos(k_b_st9, n + 1, n + 1, "PairOxrna2Stk:b_st9");
  memoryKK->create_kokkos(k_dtheta_st9_c, n + 1, n + 1, "PairOxrna2Stk:dtheta_st9_c");

  memoryKK->create_kokkos(k_a_st10, n + 1, n + 1, "PairOxrna2Stk:a_st10");
  memoryKK->create_kokkos(k_theta_st10_0, n + 1, n + 1, "PairOxrna2Stk:theta_st10_0");
  memoryKK->create_kokkos(k_dtheta_st10_ast, n + 1, n + 1, "PairOxrna2Stk:dtheta_st10_ast");
  memoryKK->create_kokkos(k_b_st10, n + 1, n + 1, "PairOxrna2Stk:b_st10");
  memoryKK->create_kokkos(k_dtheta_st10_c, n + 1, n + 1, "PairOxrna2Stk:dtheta_st10_c");

  memoryKK->create_kokkos(k_a_st1, n + 1, n + 1, "PairOxrna2Stk:a_st1");
  memoryKK->create_kokkos(k_cosphi_st1_ast, n + 1, n + 1, "PairOxrna2Stk:cosphi_st1_ast");
  memoryKK->create_kokkos(k_b_st1, n + 1, n + 1, "PairOxrna2Stk:b_st1");
  memoryKK->create_kokkos(k_cosphi_st1_c, n + 1, n + 1, "PairOxrna2Stk:cosphi_st1_c");
  memoryKK->create_kokkos(k_a_st2, n + 1, n + 1, "PairOxrna2Stk:a_st2");
  memoryKK->create_kokkos(k_cosphi_st2_ast, n + 1, n + 1, "PairOxrna2Stk:cosphi_st2_ast");
  memoryKK->create_kokkos(k_b_st2, n + 1, n + 1, "PairOxrna2Stk:b_st2");
  memoryKK->create_kokkos(k_cosphi_st2_c, n + 1, n + 1, "PairOxrna2Stk:cosphi_st2_c");

  d_epsilon_st = k_epsilon_st.template view<DeviceType>();
  d_a_st = k_a_st.template view<DeviceType>();
  d_cut_st_0 = k_cut_st_0.template view<DeviceType>();
  d_cut_st_c = k_cut_st_c.template view<DeviceType>();
  d_cut_st_lo = k_cut_st_lo.template view<DeviceType>();
  d_cut_st_hi = k_cut_st_hi.template view<DeviceType>();
  d_cut_st_lc = k_cut_st_lc.template view<DeviceType>();
  d_cut_st_hc = k_cut_st_hc.template view<DeviceType>();
  d_b_st_lo = k_b_st_lo.template view<DeviceType>();
  d_b_st_hi = k_b_st_hi.template view<DeviceType>();
  d_shift_st = k_shift_st.template view<DeviceType>();
  d_cutsq_st_hc = k_cutsq_st_hc.template view<DeviceType>();

  d_a_st5 = k_a_st5.template view<DeviceType>();
  d_theta_st5_0 = k_theta_st5_0.template view<DeviceType>();
  d_dtheta_st5_ast = k_dtheta_st5_ast.template view<DeviceType>();
  d_b_st5 = k_b_st5.template view<DeviceType>();
  d_dtheta_st5_c = k_dtheta_st5_c.template view<DeviceType>();

  d_a_st6 = k_a_st6.template view<DeviceType>();
  d_theta_st6_0 = k_theta_st6_0.template view<DeviceType>();
  d_dtheta_st6_ast = k_dtheta_st6_ast.template view<DeviceType>();
  d_b_st6 = k_b_st6.template view<DeviceType>();
  d_dtheta_st6_c = k_dtheta_st6_c.template view<DeviceType>();

  d_a_st9 = k_a_st9.template view<DeviceType>();
  d_theta_st9_0 = k_theta_st9_0.template view<DeviceType>();
  d_dtheta_st9_ast = k_dtheta_st9_ast.template view<DeviceType>();
  d_b_st9 = k_b_st9.template view<DeviceType>();
  d_dtheta_st9_c = k_dtheta_st9_c.template view<DeviceType>();

  d_a_st10 = k_a_st10.template view<DeviceType>();
  d_theta_st10_0 = k_theta_st10_0.template view<DeviceType>();
  d_dtheta_st10_ast = k_dtheta_st10_ast.template view<DeviceType>();
  d_b_st10 = k_b_st10.template view<DeviceType>();
  d_dtheta_st10_c = k_dtheta_st10_c.template view<DeviceType>();

  d_a_st1 = k_a_st1.template view<DeviceType>();
  d_cosphi_st1_ast = k_cosphi_st1_ast.template view<DeviceType>();
  d_b_st1 = k_b_st1.template view<DeviceType>();
  d_cosphi_st1_c = k_cosphi_st1_c.template view<DeviceType>();
  d_a_st2 = k_a_st2.template view<DeviceType>();
  d_cosphi_st2_ast = k_cosphi_st2_ast.template view<DeviceType>();
  d_b_st2 = k_b_st2.template view<DeviceType>();
  d_cosphi_st2_c = k_cosphi_st2_c.template view<DeviceType>();
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2StkKokkos<DeviceType>::settings(int narg, char ** /*arg*/)
{
  if (narg != 0) error->all(FLERR, "Illegal pair_style command");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2StkKokkos<DeviceType>::init_style()
{
  if (!atom->style_match("oxdna")) {
    error->all(FLERR,
               "Must use 'atom_style hybrid bond ellipsoid oxdna' with pair style "
               "oxrna2/stk/kk");
  }

  neighbor->add_request(this);
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType, LMPHostType> &&
                           !std::is_same_v<DeviceType, LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType, LMPDeviceType>);

  fix_oxdna_lrfKK = nullptr;
  auto fixes = modify->get_fix_by_style("^OXDNA/LRF/kk");
  if (fixes.size() == 0)
    error->all(FLERR, "Fix OXDNA/LRF/kk not found. Ensure pair ox*na*/excv/kk is present");
  else
    fix_oxdna_lrfKK = dynamic_cast<FixOxdnaLRFKokkos<DeviceType> *>(fixes[0]);

  auto prime_fixes = modify->get_fix_by_style("^OXDNA/PRIME_NEIGHS/kk");
  if (prime_fixes.size() == 0)
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(
        modify->add_fix("prime_neighs_kk all OXDNA/PRIME_NEIGHS/kk"));
  else
    fix_oxdna_prime_neighsKK = dynamic_cast<FixOxdnaPrimeNeighsKokkos<DeviceType> *>(prime_fixes[0]);

  if (!fix_oxdna_prime_neighsKK) error->all(FLERR, "Fix OXDNA/PRIME_NEIGHS/kk not found");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
double PairOxrna2StkKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairOxrna2Stk::init_one(i, j);

  k_epsilon_st.view_host()(i, j) = epsilon_st[i][j];
  k_epsilon_st.view_host()(j, i) = epsilon_st[j][i];
  k_a_st.view_host()(i, j) = a_st[i][j];
  k_a_st.view_host()(j, i) = a_st[j][i];
  k_cut_st_0.view_host()(i, j) = cut_st_0[i][j];
  k_cut_st_0.view_host()(j, i) = cut_st_0[j][i];
  k_cut_st_c.view_host()(i, j) = cut_st_c[i][j];
  k_cut_st_c.view_host()(j, i) = cut_st_c[j][i];
  k_cut_st_lo.view_host()(i, j) = cut_st_lo[i][j];
  k_cut_st_lo.view_host()(j, i) = cut_st_lo[j][i];
  k_cut_st_hi.view_host()(i, j) = cut_st_hi[i][j];
  k_cut_st_hi.view_host()(j, i) = cut_st_hi[j][i];
  k_cut_st_lc.view_host()(i, j) = cut_st_lc[i][j];
  k_cut_st_lc.view_host()(j, i) = cut_st_lc[j][i];
  k_cut_st_hc.view_host()(i, j) = cut_st_hc[i][j];
  k_cut_st_hc.view_host()(j, i) = cut_st_hc[j][i];
  k_b_st_lo.view_host()(i, j) = b_st_lo[i][j];
  k_b_st_lo.view_host()(j, i) = b_st_lo[j][i];
  k_b_st_hi.view_host()(i, j) = b_st_hi[i][j];
  k_b_st_hi.view_host()(j, i) = b_st_hi[j][i];
  k_shift_st.view_host()(i, j) = shift_st[i][j];
  k_shift_st.view_host()(j, i) = shift_st[j][i];
  k_cutsq_st_hc.view_host()(i, j) = cutsq_st_hc[i][j];
  k_cutsq_st_hc.view_host()(j, i) = cutsq_st_hc[j][i];

  k_a_st5.view_host()(i, j) = a_st5[i][j];
  k_a_st5.view_host()(j, i) = a_st5[j][i];
  k_theta_st5_0.view_host()(i, j) = theta_st5_0[i][j];
  k_theta_st5_0.view_host()(j, i) = theta_st5_0[j][i];
  k_dtheta_st5_ast.view_host()(i, j) = dtheta_st5_ast[i][j];
  k_dtheta_st5_ast.view_host()(j, i) = dtheta_st5_ast[j][i];
  k_b_st5.view_host()(i, j) = b_st5[i][j];
  k_b_st5.view_host()(j, i) = b_st5[j][i];
  k_dtheta_st5_c.view_host()(i, j) = dtheta_st5_c[i][j];
  k_dtheta_st5_c.view_host()(j, i) = dtheta_st5_c[j][i];

  k_a_st6.view_host()(i, j) = a_st6[i][j];
  k_a_st6.view_host()(j, i) = a_st6[j][i];
  k_theta_st6_0.view_host()(i, j) = theta_st6_0[i][j];
  k_theta_st6_0.view_host()(j, i) = theta_st6_0[j][i];
  k_dtheta_st6_ast.view_host()(i, j) = dtheta_st6_ast[i][j];
  k_dtheta_st6_ast.view_host()(j, i) = dtheta_st6_ast[j][i];
  k_b_st6.view_host()(i, j) = b_st6[i][j];
  k_b_st6.view_host()(j, i) = b_st6[j][i];
  k_dtheta_st6_c.view_host()(i, j) = dtheta_st6_c[i][j];
  k_dtheta_st6_c.view_host()(j, i) = dtheta_st6_c[j][i];

  k_a_st9.view_host()(i, j) = a_st9[i][j];
  k_a_st9.view_host()(j, i) = a_st9[j][i];
  k_theta_st9_0.view_host()(i, j) = theta_st9_0[i][j];
  k_theta_st9_0.view_host()(j, i) = theta_st9_0[j][i];
  k_dtheta_st9_ast.view_host()(i, j) = dtheta_st9_ast[i][j];
  k_dtheta_st9_ast.view_host()(j, i) = dtheta_st9_ast[j][i];
  k_b_st9.view_host()(i, j) = b_st9[i][j];
  k_b_st9.view_host()(j, i) = b_st9[j][i];
  k_dtheta_st9_c.view_host()(i, j) = dtheta_st9_c[i][j];
  k_dtheta_st9_c.view_host()(j, i) = dtheta_st9_c[j][i];

  k_a_st10.view_host()(i, j) = a_st10[i][j];
  k_a_st10.view_host()(j, i) = a_st10[j][i];
  k_theta_st10_0.view_host()(i, j) = theta_st10_0[i][j];
  k_theta_st10_0.view_host()(j, i) = theta_st10_0[j][i];
  k_dtheta_st10_ast.view_host()(i, j) = dtheta_st10_ast[i][j];
  k_dtheta_st10_ast.view_host()(j, i) = dtheta_st10_ast[j][i];
  k_b_st10.view_host()(i, j) = b_st10[i][j];
  k_b_st10.view_host()(j, i) = b_st10[j][i];
  k_dtheta_st10_c.view_host()(i, j) = dtheta_st10_c[i][j];
  k_dtheta_st10_c.view_host()(j, i) = dtheta_st10_c[j][i];

  k_a_st1.view_host()(i, j) = a_st1[i][j];
  k_a_st1.view_host()(j, i) = a_st1[j][i];
  k_cosphi_st1_ast.view_host()(i, j) = cosphi_st1_ast[i][j];
  k_cosphi_st1_ast.view_host()(j, i) = cosphi_st1_ast[j][i];
  k_b_st1.view_host()(i, j) = b_st1[i][j];
  k_b_st1.view_host()(j, i) = b_st1[j][i];
  k_cosphi_st1_c.view_host()(i, j) = cosphi_st1_c[i][j];
  k_cosphi_st1_c.view_host()(j, i) = cosphi_st1_c[j][i];

  k_a_st2.view_host()(i, j) = a_st2[i][j];
  k_a_st2.view_host()(j, i) = a_st2[j][i];
  k_cosphi_st2_ast.view_host()(i, j) = cosphi_st2_ast[i][j];
  k_cosphi_st2_ast.view_host()(j, i) = cosphi_st2_ast[j][i];
  k_b_st2.view_host()(i, j) = b_st2[i][j];
  k_b_st2.view_host()(j, i) = b_st2[j][i];
  k_cosphi_st2_c.view_host()(i, j) = cosphi_st2_c[i][j];
  k_cosphi_st2_c.view_host()(j, i) = cosphi_st2_c[j][i];

  k_epsilon_st.template modify<LMPHostType>();
  k_a_st.template modify<LMPHostType>();
  k_cut_st_0.template modify<LMPHostType>();
  k_cut_st_c.template modify<LMPHostType>();
  k_cut_st_lo.template modify<LMPHostType>();
  k_cut_st_hi.template modify<LMPHostType>();
  k_cut_st_lc.template modify<LMPHostType>();
  k_cut_st_hc.template modify<LMPHostType>();
  k_b_st_lo.template modify<LMPHostType>();
  k_b_st_hi.template modify<LMPHostType>();
  k_shift_st.template modify<LMPHostType>();
  k_cutsq_st_hc.template modify<LMPHostType>();

  k_a_st5.template modify<LMPHostType>();
  k_theta_st5_0.template modify<LMPHostType>();
  k_dtheta_st5_ast.template modify<LMPHostType>();
  k_b_st5.template modify<LMPHostType>();
  k_dtheta_st5_c.template modify<LMPHostType>();

  k_a_st6.template modify<LMPHostType>();
  k_theta_st6_0.template modify<LMPHostType>();
  k_dtheta_st6_ast.template modify<LMPHostType>();
  k_b_st6.template modify<LMPHostType>();
  k_dtheta_st6_c.template modify<LMPHostType>();

  k_a_st9.template modify<LMPHostType>();
  k_theta_st9_0.template modify<LMPHostType>();
  k_dtheta_st9_ast.template modify<LMPHostType>();
  k_b_st9.template modify<LMPHostType>();
  k_dtheta_st9_c.template modify<LMPHostType>();

  k_a_st10.template modify<LMPHostType>();
  k_theta_st10_0.template modify<LMPHostType>();
  k_dtheta_st10_ast.template modify<LMPHostType>();
  k_b_st10.template modify<LMPHostType>();
  k_dtheta_st10_c.template modify<LMPHostType>();

  k_a_st1.template modify<LMPHostType>();
  k_cosphi_st1_ast.template modify<LMPHostType>();
  k_b_st1.template modify<LMPHostType>();
  k_cosphi_st1_c.template modify<LMPHostType>();

  k_a_st2.template modify<LMPHostType>();
  k_cosphi_st2_ast.template modify<LMPHostType>();
  k_b_st2.template modify<LMPHostType>();
  k_cosphi_st2_c.template modify<LMPHostType>();

  k_epsilon_st.template sync<DeviceType>();
  k_a_st.template sync<DeviceType>();
  k_cut_st_0.template sync<DeviceType>();
  k_cut_st_c.template sync<DeviceType>();
  k_cut_st_lo.template sync<DeviceType>();
  k_cut_st_hi.template sync<DeviceType>();
  k_cut_st_lc.template sync<DeviceType>();
  k_cut_st_hc.template sync<DeviceType>();
  k_b_st_lo.template sync<DeviceType>();
  k_b_st_hi.template sync<DeviceType>();
  k_shift_st.template sync<DeviceType>();
  k_cutsq_st_hc.template sync<DeviceType>();

  k_a_st5.template sync<DeviceType>();
  k_theta_st5_0.template sync<DeviceType>();
  k_dtheta_st5_ast.template sync<DeviceType>();
  k_b_st5.template sync<DeviceType>();
  k_dtheta_st5_c.template sync<DeviceType>();

  k_a_st6.template sync<DeviceType>();
  k_theta_st6_0.template sync<DeviceType>();
  k_dtheta_st6_ast.template sync<DeviceType>();
  k_b_st6.template sync<DeviceType>();
  k_dtheta_st6_c.template sync<DeviceType>();

  k_a_st9.template sync<DeviceType>();
  k_theta_st9_0.template sync<DeviceType>();
  k_dtheta_st9_ast.template sync<DeviceType>();
  k_b_st9.template sync<DeviceType>();
  k_dtheta_st9_c.template sync<DeviceType>();

  k_a_st10.template sync<DeviceType>();
  k_theta_st10_0.template sync<DeviceType>();
  k_dtheta_st10_ast.template sync<DeviceType>();
  k_b_st10.template sync<DeviceType>();
  k_dtheta_st10_c.template sync<DeviceType>();

  k_a_st1.template sync<DeviceType>();
  k_cosphi_st1_ast.template sync<DeviceType>();
  k_b_st1.template sync<DeviceType>();
  k_cosphi_st1_c.template sync<DeviceType>();

  k_a_st2.template sync<DeviceType>();
  k_cosphi_st2_ast.template sync<DeviceType>();
  k_b_st2.template sync<DeviceType>();
  k_cosphi_st2_c.template sync<DeviceType>();

  return cutone;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairOxrna2StkKokkos<DeviceType>::coeff(int narg, char **arg)
{
  PairOxrna2Stk::coeff(narg, arg);
}

/* ----------------------------------------------------------------------
   tally energy and virial into global and per-atom accumulators

   NOTE: Although this is a pair style interaction, the algorithm below
   follows the virial incrementation of the bond style. This is because
   the bond topology is used in the main compute loop.
------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION void PairOxrna2StkKokkos<DeviceType>::ev_tally_xyz(
    EV_FLOAT &ev, const int &i, const int &j, const int &nlocal, const int &newton_bond,
    const KK_FLOAT &evdwl, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy,
    const KK_ACC_FLOAT &fz, const KK_FLOAT &delx, const KK_FLOAT &dely,
    const KK_FLOAT &delz) const
{
  KK_ACC_FLOAT evdwlhalf;
  KK_ACC_FLOAT v[6];

  // The eatom and vatom arrays are atomic.
  Kokkos::View<KK_ACC_FLOAT *, typename DAT::t_kkacc_1d::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      v_eatom = d_eatom;
  Kokkos::View<KK_ACC_FLOAT *[6], typename DAT::t_kkacc_1d_6::array_layout,
               typename KKDevice<DeviceType>::value,
               Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Unmanaged>>
      v_vatom = d_vatom;

  if (eflag_either) {
    if (eflag_global) {
      if (newton_bond)
        ev.evdwl += static_cast<KK_ACC_FLOAT>(evdwl);
      else {
        evdwlhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * evdwl);
        if (i < nlocal) ev.evdwl += evdwlhalf;
        if (j < nlocal) ev.evdwl += evdwlhalf;
      }
    }
    if (eflag_atom) {
      evdwlhalf = static_cast<KK_ACC_FLOAT>(static_cast<KK_FLOAT>(0.5) * evdwl);
      if (newton_bond || i < nlocal) v_eatom[i] += evdwlhalf;
      if (newton_bond || j < nlocal) v_eatom[j] += evdwlhalf;
    }
  }

  if (vflag_either) {
    v[0] = static_cast<KK_ACC_FLOAT>(delx * fx);
    v[1] = static_cast<KK_ACC_FLOAT>(dely * fy);
    v[2] = static_cast<KK_ACC_FLOAT>(delz * fz);
    v[3] = static_cast<KK_ACC_FLOAT>(delx * fy);
    v[4] = static_cast<KK_ACC_FLOAT>(delx * fz);
    v[5] = static_cast<KK_ACC_FLOAT>(dely * fz);

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
          ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5 * v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5 * v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5 * v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5 * v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5 * v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5 * v[5]);
        }
        if (j < nlocal) {
          ev.v[0] += static_cast<KK_ACC_FLOAT>(0.5 * v[0]);
          ev.v[1] += static_cast<KK_ACC_FLOAT>(0.5 * v[1]);
          ev.v[2] += static_cast<KK_ACC_FLOAT>(0.5 * v[2]);
          ev.v[3] += static_cast<KK_ACC_FLOAT>(0.5 * v[3]);
          ev.v[4] += static_cast<KK_ACC_FLOAT>(0.5 * v[4]);
          ev.v[5] += static_cast<KK_ACC_FLOAT>(0.5 * v[5]);
        }
      }
    }

    if (vflag_atom) {
      if (newton_bond || i < nlocal) {
        v_vatom(i, 0) += static_cast<KK_ACC_FLOAT>(0.5 * v[0]);
        v_vatom(i, 1) += static_cast<KK_ACC_FLOAT>(0.5 * v[1]);
        v_vatom(i, 2) += static_cast<KK_ACC_FLOAT>(0.5 * v[2]);
        v_vatom(i, 3) += static_cast<KK_ACC_FLOAT>(0.5 * v[3]);
        v_vatom(i, 4) += static_cast<KK_ACC_FLOAT>(0.5 * v[4]);
        v_vatom(i, 5) += static_cast<KK_ACC_FLOAT>(0.5 * v[5]);
      }
      if (newton_bond || j < nlocal) {
        v_vatom(j, 0) += static_cast<KK_ACC_FLOAT>(0.5 * v[0]);
        v_vatom(j, 1) += static_cast<KK_ACC_FLOAT>(0.5 * v[1]);
        v_vatom(j, 2) += static_cast<KK_ACC_FLOAT>(0.5 * v[2]);
        v_vatom(j, 3) += static_cast<KK_ACC_FLOAT>(0.5 * v[3]);
        v_vatom(j, 4) += static_cast<KK_ACC_FLOAT>(0.5 * v[4]);
        v_vatom(j, 5) += static_cast<KK_ACC_FLOAT>(0.5 * v[5]);
      }
    }
  }
}

namespace LAMMPS_NS {
template class PairOxrna2StkKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairOxrna2StkKokkos<LMPHostType>;
#endif
}    // namespace LAMMPS_NS
