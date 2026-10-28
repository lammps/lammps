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

// Device functions of the screened-pair kernel of PairOxdnaHbondKokkos.
// They are in a header so that the fused hbond + oxdna3/xstk kernel in
// pair_oxdna3_xstk_kokkos.cpp can inline them (GPU builds do not use
// relocatable device code).

#ifndef LMP_PAIR_OXDNA_HBOND_KOKKOS_IMPL_H
#define LMP_PAIR_OXDNA_HBOND_KOKKOS_IMPL_H

#include "pair_oxdna_hbond_kokkos.h"

#include "mf_oxdna_kokkos.h"

namespace LAMMPS_NS {

using namespace MFOxdnaKokkos;

/* ----------------------------------------------------------------------
   ComputeGPUPair Functor(s) and staged hbond helpers for lower
   live register pressure in GPU kernels.
-------------------------------------------------------------------------- */

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_radial_terms(const int &atype, const int &btype,
  const KK_FLOAT &r_hb, KK_FLOAT &f1, KK_FLOAT &df1) const
{
  const KK_FLOAT p_epsilon_hb = d_params_hb(atype,btype).epsilon_hb;
  const KK_FLOAT p_a_hb = d_params_hb(atype,btype).a_hb;
  const KK_FLOAT p_cut_hb_0 = d_params_hb(atype,btype).cut_hb_0;
  const KK_FLOAT p_cut_hb_lc = d_params_hb(atype,btype).cut_hb_lc;
  const KK_FLOAT p_cut_hb_hc = d_params_hb(atype,btype).cut_hb_hc;
  const KK_FLOAT p_cut_hb_lo = d_params_hb(atype,btype).cut_hb_lo;
  const KK_FLOAT p_cut_hb_hi = d_params_hb(atype,btype).cut_hb_hi;
  const KK_FLOAT p_b_hb_lo = d_params_hb(atype,btype).b_hb_lo;
  const KK_FLOAT p_b_hb_hi = d_params_hb(atype,btype).b_hb_hi;
  const KK_FLOAT p_shift_hb = d_params_hb(atype,btype).shift_hb;

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
  const KK_FLOAT p_a_hb1 = d_params_hb(atype,btype).a_hb1;
  const KK_FLOAT p_theta_hb1_0 = d_params_hb(atype,btype).theta_hb1_0;
  const KK_FLOAT p_dtheta_hb1_ast = d_params_hb(atype,btype).dtheta_hb1_ast;
  const KK_FLOAT p_b_hb1 = d_params_hb(atype,btype).b_hb1;
  const KK_FLOAT p_dtheta_hb1_c = d_params_hb(atype,btype).dtheta_hb1_c;

  KK_FLOAT cost1 = -Kokkos::fma(a_nx[2], b_nx[2], Kokkos::fma(a_nx[1], b_nx[1], a_nx[0] * b_nx[0]));
  if (cost1 > static_cast<KK_FLOAT>(1.0)) cost1 = static_cast<KK_FLOAT>(1.0);
  if (cost1 < static_cast<KK_FLOAT>(-1.0)) cost1 = static_cast<KK_FLOAT>(-1.0);
  theta1 = Kokkos::acos(cost1);

  f4t1 = F4_KK(theta1, p_a_hb1, p_theta_hb1_0, p_dtheta_hb1_ast, p_b_hb1, p_dtheta_hb1_c);
  if (f4t1 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin1_sq = Kokkos::fma(-cost1, cost1, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin1 = (sin1_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin1_sq) : static_cast<KK_FLOAT>(0.0);
  df4t1 = DF4_KK(theta1, p_a_hb1, p_theta_hb1_0, p_dtheta_hb1_ast, p_b_hb1, p_dtheta_hb1_c) * rsin1;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta2_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta2, KK_FLOAT &cost2, KK_FLOAT &f4t2, KK_FLOAT &df4t2) const
{
  const KK_FLOAT p_a_hb2 = d_params_hb(atype,btype).a_hb2;
  const KK_FLOAT p_theta_hb2_0 = d_params_hb(atype,btype).theta_hb2_0;
  const KK_FLOAT p_dtheta_hb2_ast = d_params_hb(atype,btype).dtheta_hb2_ast;
  const KK_FLOAT p_b_hb2 = d_params_hb(atype,btype).b_hb2;
  const KK_FLOAT p_dtheta_hb2_c = d_params_hb(atype,btype).dtheta_hb2_c;

  cost2 = -Kokkos::fma(a_nx[2], delr_hb_norm[2], Kokkos::fma(a_nx[1], delr_hb_norm[1], a_nx[0] * delr_hb_norm[0]));
  if (cost2 > static_cast<KK_FLOAT>(1.0)) cost2 = static_cast<KK_FLOAT>(1.0);
  if (cost2 < static_cast<KK_FLOAT>(-1.0)) cost2 = static_cast<KK_FLOAT>(-1.0);
  theta2 = Kokkos::acos(cost2);

  f4t2 = F4_KK(theta2, p_a_hb2, p_theta_hb2_0, p_dtheta_hb2_ast, p_b_hb2, p_dtheta_hb2_c);
  if (f4t2 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin2_sq = Kokkos::fma(-cost2, cost2, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin2 = (sin2_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin2_sq) : static_cast<KK_FLOAT>(0.0);
  df4t2 = DF4_KK(theta2, p_a_hb2, p_theta_hb2_0, p_dtheta_hb2_ast, p_b_hb2, p_dtheta_hb2_c) * rsin2;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta3_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nx)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta3, KK_FLOAT &cost3, KK_FLOAT &f4t3, KK_FLOAT &df4t3) const
{
  const KK_FLOAT p_a_hb3 = d_params_hb(atype,btype).a_hb3;
  const KK_FLOAT p_theta_hb3_0 = d_params_hb(atype,btype).theta_hb3_0;
  const KK_FLOAT p_dtheta_hb3_ast = d_params_hb(atype,btype).dtheta_hb3_ast;
  const KK_FLOAT p_b_hb3 = d_params_hb(atype,btype).b_hb3;
  const KK_FLOAT p_dtheta_hb3_c = d_params_hb(atype,btype).dtheta_hb3_c;

  cost3 = Kokkos::fma(b_nx[2], delr_hb_norm[2], Kokkos::fma(b_nx[1], delr_hb_norm[1], b_nx[0] * delr_hb_norm[0]));
  if (cost3 > static_cast<KK_FLOAT>(1.0)) cost3 = static_cast<KK_FLOAT>(1.0);
  if (cost3 < static_cast<KK_FLOAT>(-1.0)) cost3 = static_cast<KK_FLOAT>(-1.0);
  theta3 = Kokkos::acos(cost3);

  f4t3 = F4_KK(theta3, p_a_hb3, p_theta_hb3_0, p_dtheta_hb3_ast, p_b_hb3, p_dtheta_hb3_c);
  if (f4t3 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin3_sq = Kokkos::fma(-cost3, cost3, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin3 = (sin3_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin3_sq) : static_cast<KK_FLOAT>(0.0);
  df4t3 = DF4_KK(theta3, p_a_hb3, p_theta_hb3_0, p_dtheta_hb3_ast, p_b_hb3, p_dtheta_hb3_c) * rsin3;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta4_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&b_nz)[3],
  KK_FLOAT &theta4, KK_FLOAT &f4t4, KK_FLOAT &df4t4) const
{
  const KK_FLOAT p_a_hb4 = d_params_hb(atype,btype).a_hb4;
  const KK_FLOAT p_theta_hb4_0 = d_params_hb(atype,btype).theta_hb4_0;
  const KK_FLOAT p_dtheta_hb4_ast = d_params_hb(atype,btype).dtheta_hb4_ast;
  const KK_FLOAT p_b_hb4 = d_params_hb(atype,btype).b_hb4;
  const KK_FLOAT p_dtheta_hb4_c = d_params_hb(atype,btype).dtheta_hb4_c;

  KK_FLOAT cost4 = Kokkos::fma(a_nz[2], b_nz[2], Kokkos::fma(a_nz[1], b_nz[1], a_nz[0] * b_nz[0]));
  if (cost4 > static_cast<KK_FLOAT>(1.0)) cost4 = static_cast<KK_FLOAT>(1.0);
  if (cost4 < static_cast<KK_FLOAT>(-1.0)) cost4 = static_cast<KK_FLOAT>(-1.0);
  theta4 = Kokkos::acos(cost4);

  f4t4 = F4_KK(theta4, p_a_hb4, p_theta_hb4_0, p_dtheta_hb4_ast, p_b_hb4, p_dtheta_hb4_c);
  if (f4t4 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin4_sq = Kokkos::fma(-cost4, cost4, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin4 = (sin4_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin4_sq) : static_cast<KK_FLOAT>(0.0);
  df4t4 = DF4_KK(theta4, p_a_hb4, p_theta_hb4_0, p_dtheta_hb4_ast, p_b_hb4, p_dtheta_hb4_c) * rsin4;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta7_terms(const int &atype, const int &btype,
  const KK_FLOAT (&a_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta7, KK_FLOAT &cost7, KK_FLOAT &f4t7, KK_FLOAT &df4t7) const
{
  const KK_FLOAT p_a_hb7 = d_params_hb(atype,btype).a_hb7;
  const KK_FLOAT p_theta_hb7_0 = d_params_hb(atype,btype).theta_hb7_0;
  const KK_FLOAT p_dtheta_hb7_ast = d_params_hb(atype,btype).dtheta_hb7_ast;
  const KK_FLOAT p_b_hb7 = d_params_hb(atype,btype).b_hb7;
  const KK_FLOAT p_dtheta_hb7_c = d_params_hb(atype,btype).dtheta_hb7_c;

  cost7 = -Kokkos::fma(a_nz[2], delr_hb_norm[2], Kokkos::fma(a_nz[1], delr_hb_norm[1], a_nz[0] * delr_hb_norm[0]));
  if (cost7 > static_cast<KK_FLOAT>(1.0)) cost7 = static_cast<KK_FLOAT>(1.0);
  if (cost7 < static_cast<KK_FLOAT>(-1.0)) cost7 = static_cast<KK_FLOAT>(-1.0);
  theta7 = Kokkos::acos(cost7);

  f4t7 = F4_KK(theta7, p_a_hb7, p_theta_hb7_0, p_dtheta_hb7_ast, p_b_hb7, p_dtheta_hb7_c);
  if (f4t7 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin7_sq = Kokkos::fma(-cost7, cost7, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin7 = (sin7_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin7_sq) : static_cast<KK_FLOAT>(0.0);
  df4t7 = DF4_KK(theta7, p_a_hb7, p_theta_hb7_0, p_dtheta_hb7_ast, p_b_hb7, p_dtheta_hb7_c) * rsin7;
  return true;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::hbond_theta8_terms(const int &atype, const int &btype,
  const KK_FLOAT (&b_nz)[3], const KK_FLOAT (&delr_hb_norm)[3],
  KK_FLOAT &theta8, KK_FLOAT &cost8, KK_FLOAT &f4t8, KK_FLOAT &df4t8) const
{
  const KK_FLOAT p_a_hb8 = d_params_hb(atype,btype).a_hb8;
  const KK_FLOAT p_theta_hb8_0 = d_params_hb(atype,btype).theta_hb8_0;
  const KK_FLOAT p_dtheta_hb8_ast = d_params_hb(atype,btype).dtheta_hb8_ast;
  const KK_FLOAT p_b_hb8 = d_params_hb(atype,btype).b_hb8;
  const KK_FLOAT p_dtheta_hb8_c = d_params_hb(atype,btype).dtheta_hb8_c;

  cost8 = Kokkos::fma(b_nz[2], delr_hb_norm[2], Kokkos::fma(b_nz[1], delr_hb_norm[1], b_nz[0] * delr_hb_norm[0]));
  if (cost8 > static_cast<KK_FLOAT>(1.0)) cost8 = static_cast<KK_FLOAT>(1.0);
  if (cost8 < static_cast<KK_FLOAT>(-1.0)) cost8 = static_cast<KK_FLOAT>(-1.0);
  theta8 = Kokkos::acos(cost8);

  f4t8 = F4_KK(theta8, p_a_hb8, p_theta_hb8_0, p_dtheta_hb8_ast, p_b_hb8, p_dtheta_hb8_c);
  if (f4t8 == static_cast<KK_FLOAT>(0.0)) return false;

  KK_FLOAT sin8_sq = Kokkos::fma(-cost8, cost8, static_cast<KK_FLOAT>(1.0));
  // at sin(theta) = 0 the angular force and torque directions vanish, so the
  // derivative term is zero, but the pair still contributes its energy and
  // its other force terms
  const KK_FLOAT rsin8 = (sin8_sq > static_cast<KK_FLOAT>(0.0)) ?
    static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(sin8_sq) : static_cast<KK_FLOAT>(0.0);
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

  delta[0] = Kokkos::fma(ra_chb[1], delf[2], -static_cast<KK_ACC_FLOAT>(ra_chb[2]) * delf[1]);
  delta[1] = Kokkos::fma(ra_chb[2], delf[0], -static_cast<KK_ACC_FLOAT>(ra_chb[0]) * delf[2]);
  delta[2] = Kokkos::fma(ra_chb[0], delf[1], -static_cast<KK_ACC_FLOAT>(ra_chb[1]) * delf[0]);

  deltb[0] = Kokkos::fma(rb_chb[1], delf[2], -static_cast<KK_ACC_FLOAT>(rb_chb[2]) * delf[1]);
  deltb[1] = Kokkos::fma(rb_chb[2], delf[0], -static_cast<KK_ACC_FLOAT>(rb_chb[0]) * delf[2]);
  deltb[2] = Kokkos::fma(rb_chb[0], delf[1], -static_cast<KK_ACC_FLOAT>(rb_chb[1]) * delf[0]);
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
template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG, int RADIAL_ONLY>
KOKKOS_INLINE_FUNCTION
bool PairOxdnaHbondKokkos<DeviceType>::screened_pair_body(TagPairOxdnaHbondComputeGPUPair<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int &ipair,
  KK_ACC_FLOAT (&fa)[3], KK_ACC_FLOAT (&ta)[3], EV_FLOAT &ev) const
{
  KK_ACC_FLOAT rxf_a[3], rxf_b[3] = {0.0, 0.0, 0.0};    // r x f torques on a and b
  // one thread per neighbor pair: several threads update the same atoms
  // with any neighbor list style, so all updates must be atomic

  auto v_f = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,decltype(dup_f),decltype(ndup_f)>::get(dup_f,ndup_f);
  auto a_f = v_f.template access<Kokkos::Experimental::ScatterAtomic>();
  auto v_torque = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_torque),decltype(ndup_torque)>::get(dup_torque,ndup_torque);
  auto a_torque = v_torque.template access<Kokkos::Experimental::ScatterAtomic>();

  // Direct packed pair lookup: high 32 bits = a, low 32 bits = b.
  const uint64_t pair = d_pairs_screened(ipair);
  // "pair >> 32" shifts the pair to the right by 32 bits, so the upper 32 bits
  // becomes the lower 32 bits to recover the atom-a index.
  const int a = static_cast<int>(pair >> 32);
  const int atype = type(a);
  // "pair & 0xffffffffu" keeps only the lower 32 bits to recover the atom-b index.
  int b = static_cast<int>(pair & 0xffffffffu);
  const KK_FLOAT factor_lj = static_cast<KK_FLOAT>(special_lj[sbmask(b)]);
  if (factor_lj == static_cast<KK_FLOAT>(0.0)) return false;
  b &= NEIGHMASK;
  const int btype = type(b);

  // no hydrogen bonding between these base types (e.g. non-complementary bases):
  // f1 and thus the energy would be zero, so skip the site geometry altogether
  if (d_params_hb(atype,btype).epsilon_hb == static_cast<KK_FLOAT>(0.0)) return false;

  if (unique_basepair_enabled) {
    const int idca = d_idc(a);
    const int idcb = d_idc(b);
    if (idca != tag(b) && idcb != tag(a) && idca > 0 && idcb > 0) return false;
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
    constexpr KK_FLOAT dx_cbs_oxdna1=static_cast<KK_FLOAT>(+0.4);
    ra_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(a,0);
    ra_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(a,1);
    ra_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(a,2);
    rb_chb[0] = dx_cbs_oxdna1*d_nx_xtrct(b,0);
    rb_chb[1] = dx_cbs_oxdna1*d_nx_xtrct(b,1);
    rb_chb[2] = dx_cbs_oxdna1*d_nx_xtrct(b,2);
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
  if (rsq_hb <= static_cast<KK_FLOAT>(0.0)) return false;
  rinv_hb = static_cast<KK_FLOAT>(1.0) / Kokkos::sqrt(rsq_hb);
  r_hb = rsq_hb * rinv_hb;

  delr_hb_norm[0] = delr_hb[0] * rinv_hb;
  delr_hb_norm[1] = delr_hb[1] * rinv_hb;
  delr_hb_norm[2] = delr_hb[2] * rinv_hb;

  KK_FLOAT f1, f4t1, f4t2, f4t3, f4t4, f4t7, f4t8;
  KK_FLOAT df1, df4t1, df4t2, df4t3, df4t4, df4t7, df4t8;
  KK_FLOAT theta1, theta2, theta3, theta4, theta7, theta8;
  KK_FLOAT cost2, cost3, cost7, cost8;

  if (!hbond_radial_terms(atype, btype, r_hb, f1, df1)) return false;
  if constexpr (RADIAL_ONLY) return true;
  if (!hbond_theta1_terms(atype, btype, a_nx, b_nx, theta1, f4t1, df4t1)) return false;
  if (!hbond_theta2_terms(atype, btype, a_nx, delr_hb_norm, theta2, cost2, f4t2, df4t2)) return false;
  if (!hbond_theta3_terms(atype, btype, b_nx, delr_hb_norm, theta3, cost3, f4t3, df4t3)) return false;
  if (!hbond_theta4_terms(atype, btype, a_nz, b_nz, theta4, f4t4, df4t4)) return false;
  if (!hbond_theta7_terms(atype, btype, a_nz, delr_hb_norm, theta7, cost7, f4t7, df4t7)) return false;
  if (!hbond_theta8_terms(atype, btype, b_nz, delr_hb_norm, theta8, cost8, f4t8, df4t8)) return false;

  evdwl = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8 * factor_lj;
  if (evdwl == static_cast<KK_FLOAT>(0.0)) return false;

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

  fa[0] += delf[0];
  fa[1] += delf[1];
  fa[2] += delf[2];
  // keep the r x f torques; applied together with the pure torques below
  rxf_a[0] = delta[0];
  rxf_a[1] = delta[1];
  rxf_a[2] = delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_f(b,0) -= delf[0];
    a_f(b,1) -= delf[1];
    a_f(b,2) -= delf[2];
    rxf_b[0] = deltb[0];
    rxf_b[1] = deltb[1];
    rxf_b[2] = deltb[2];
  }

  if (EVFLAG) {
    ev.evdwl += static_cast<KK_ACC_FLOAT>((((NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD)&&(NEWTON_PAIR||(b<nlocal)))?static_cast<KK_FLOAT>(1.0):static_cast<KK_FLOAT>(0.5))*evdwl);

    if (vflag_either || eflag_atom) {
      this->template ev_tally_xyz<NEIGHFLAG,NEWTON_PAIR,1>(ev,a,b,evdwl,\
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

  ta[0] += rxf_a[0] + delta[0];
  ta[1] += rxf_a[1] + delta[1];
  ta[2] += rxf_a[2] + delta[2];

  if ( (NEIGHFLAG==HALF || NEIGHFLAG==HALFTHREAD) && (NEWTON_PAIR || b < nlocal) ) {
    a_torque(b,0) -= rxf_b[0] + deltb[0];
    a_torque(b,1) -= rxf_b[1] + deltb[1];
    a_torque(b,2) -= rxf_b[2] + deltb[2];
  }
  return true;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
template<int NEIGHFLAG, int NEWTON_PAIR, int PAIRWISE>
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
  auto a_eatom = v_eatom.template access<std::conditional_t<PAIRWISE,Kokkos::Experimental::ScatterAtomic,AtomicDup_v<NEIGHFLAG,DeviceType>>>();

  auto v_vatom = ScatterViewHelper<NeedDup_v<NEIGHFLAG,DeviceType>,\
    decltype(dup_vatom),decltype(ndup_vatom)>::get(dup_vatom,ndup_vatom);
  auto a_vatom = v_vatom.template access<std::conditional_t<PAIRWISE,Kokkos::Experimental::ScatterAtomic,AtomicDup_v<NEIGHFLAG,DeviceType>>>();

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
int PairOxdnaHbondKokkos<DeviceType>::sbmask(const int& j) const {
  return j >> SBBITS & 3;
}

}    // namespace LAMMPS_NS

#endif
