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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxrna2/stk/kk,PairOxrna2StkKokkos<LMPDeviceType>);
PairStyle(oxrna2/stk/kk/device,PairOxrna2StkKokkos<LMPDeviceType>);
PairStyle(oxrna2/stk/kk/host,PairOxrna2StkKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_STK_KOKKOS_H
#define LMP_PAIR_OXRNA2_STK_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_oxrna2_stk.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration
template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos;  // forward declaration

template<int NEWTON_BOND, int EVFLAG>
struct TagPairOxrna2StkCompute {};

template<class DeviceType>
class PairOxrna2StkKokkos : public PairOxrna2Stk, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairOxrna2StkKokkos(class LAMMPS *);
  ~PairOxrna2StkKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxrna2StkCompute<NEWTON_BOND, EVFLAG>, const int &, EV_FLOAT &) const;

  template<int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxrna2StkCompute<NEWTON_BOND, EVFLAG>, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j, const int &nlocal,
                    const int &newton_bond, const KK_FLOAT &evdwl, const KK_ACC_FLOAT &fx,
                    const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz, const KK_FLOAT &delx,
                    const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 protected:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nbondlist;
  int nlocal, newton_bond, eflag, vflag;

  // stacking interaction parameters
  typename AT::tdual_kkfloat_2d k_epsilon_st, k_a_st;
  typename AT::tdual_kkfloat_2d k_cut_st_0, k_cut_st_c, k_cut_st_lo, k_cut_st_hi;
  typename AT::tdual_kkfloat_2d k_cut_st_lc, k_cut_st_hc;
  typename AT::tdual_kkfloat_2d k_b_st_lo, k_b_st_hi;
  typename AT::tdual_kkfloat_2d k_shift_st, k_cutsq_st_hc;
  typename AT::tdual_kkfloat_2d k_a_st5, k_theta_st5_0, k_dtheta_st5_ast;
  typename AT::tdual_kkfloat_2d k_b_st5, k_dtheta_st5_c;
  typename AT::tdual_kkfloat_2d k_a_st6, k_theta_st6_0, k_dtheta_st6_ast;
  typename AT::tdual_kkfloat_2d k_b_st6, k_dtheta_st6_c;
  typename AT::tdual_kkfloat_2d k_a_st9, k_theta_st9_0, k_dtheta_st9_ast;
  typename AT::tdual_kkfloat_2d k_b_st9, k_dtheta_st9_c;
  typename AT::tdual_kkfloat_2d k_a_st10, k_theta_st10_0, k_dtheta_st10_ast;
  typename AT::tdual_kkfloat_2d k_b_st10, k_dtheta_st10_c;
  typename AT::tdual_kkfloat_2d k_a_st1, k_cosphi_st1_ast, k_b_st1, k_cosphi_st1_c;
  typename AT::tdual_kkfloat_2d k_a_st2, k_cosphi_st2_ast, k_b_st2, k_cosphi_st2_c;

  typename AT::t_kkfloat_2d_randomread d_epsilon_st, d_a_st;
  typename AT::t_kkfloat_2d_randomread d_cut_st_0, d_cut_st_c, d_cut_st_lo, d_cut_st_hi;
  typename AT::t_kkfloat_2d_randomread d_cut_st_lc, d_cut_st_hc;
  typename AT::t_kkfloat_2d_randomread d_b_st_lo, d_b_st_hi;
  typename AT::t_kkfloat_2d_randomread d_shift_st, d_cutsq_st_hc;
  typename AT::t_kkfloat_2d_randomread d_a_st5, d_theta_st5_0, d_dtheta_st5_ast;
  typename AT::t_kkfloat_2d_randomread d_b_st5, d_dtheta_st5_c;
  typename AT::t_kkfloat_2d_randomread d_a_st6, d_theta_st6_0, d_dtheta_st6_ast;
  typename AT::t_kkfloat_2d_randomread d_b_st6, d_dtheta_st6_c;
  typename AT::t_kkfloat_2d_randomread d_a_st9, d_theta_st9_0, d_dtheta_st9_ast;
  typename AT::t_kkfloat_2d_randomread d_b_st9, d_dtheta_st9_c;
  typename AT::t_kkfloat_2d_randomread d_a_st10, d_theta_st10_0, d_dtheta_st10_ast;
  typename AT::t_kkfloat_2d_randomread d_b_st10, d_dtheta_st10_c;
  typename AT::t_kkfloat_2d_randomread d_a_st1, d_cosphi_st1_ast, d_b_st1, d_cosphi_st1_c;
  typename AT::t_kkfloat_2d_randomread d_a_st2, d_cosphi_st2_ast, d_b_st2, d_cosphi_st2_c;

  // per-atom arrays for local unit vectors
  typename AT::t_kkfloat_1d_3_randomread d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  void allocate() override;

  friend void pair_virial_fdotr_compute<PairOxrna2StkKokkos>(PairOxrna2StkKokkos *);

  class NeighborKokkos *neighborKK;
  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;  // ptr to OXDNA/LRF/kk fix
  FixOxdnaPrimeNeighsKokkos<DeviceType>
      *fix_oxdna_prime_neighsKK;  // ptr to OXDNA/PRIME_NEIGHS/kk fix
  bigint last_prime_neighs_bond_lastcall;

  // Precomputed atom a/b 3'/5' directionality and atom mapping of their 3' and 5' neighbors.
  // 0-3 : atom a, atom b, id3p[a], id5p[b] for each bond.
  typename AT::t_int_1d_4_randomread d_prime_neighs_bond;
};

}    // namespace LAMMPS_NS

#endif
#endif
