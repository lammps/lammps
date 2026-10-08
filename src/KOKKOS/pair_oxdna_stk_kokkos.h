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
PairStyle(oxdna/stk/kk,PairOxdnaStkKokkos<LMPDeviceType>);
PairStyle(oxdna/stk/kk/device,PairOxdnaStkKokkos<LMPDeviceType>);
PairStyle(oxdna/stk/kk/host,PairOxdnaStkKokkos<LMPHostType>);
PairStyle(oxdna2/stk/kk,PairOxdnaStkKokkos<LMPDeviceType>);
PairStyle(oxdna2/stk/kk/device,PairOxdnaStkKokkos<LMPDeviceType>);
PairStyle(oxdna2/stk/kk/host,PairOxdnaStkKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA_STK_KOKKOS_H
#define LMP_PAIR_OXDNA_STK_KOKKOS_H

#include "kokkos_base.h"
#include "fix_oxdna_lrf_kokkos.h"
#include "pair_kokkos.h"
#include "pair_oxdna_stk.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration
template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos;  // forward declaration

template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
struct TagPairOxdnaStkCompute{};

// packed per-type-pair coefficients of PairOxdnaStkKokkos
struct ParamsOxdnaStk2 {
  KK_FLOAT epsilon_st, a_st, b_st_lo, b_st_hi, theta_st4_0, a_st5;
  KK_FLOAT theta_st5_0, dtheta_st5_ast, b_st5, dtheta_st5_c, a_st6, theta_st6_0;
  KK_FLOAT dtheta_st6_ast, b_st6, dtheta_st6_c, a_st1, cosphi_st1_ast, b_st1;
  KK_FLOAT cosphi_st1_c, a_st2, cosphi_st2_ast, b_st2, cosphi_st2_c;
};

// packed per-tetramer coefficients of PairOxdnaStkKokkos
struct ParamsOxdnaStk4 {
  KK_FLOAT cut_st_0, cut_st_c, cut_st_lo, cut_st_hi, cut_st_lc, cut_st_hc;
  KK_FLOAT shift_st, cutsq_st_hc, a_st4, dtheta_st4_ast, b_st4, dtheta_st4_c;
};

template<class DeviceType>
class PairOxdnaStkKokkos : public PairOxdnaStk, public KokkosBase {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdnaStkKokkos(class LAMMPS *);
  ~PairOxdnaStkKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaStkCompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEWTON_BOND, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaStkCompute<OXDNAFLAG,NEWTON_BOND,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j, const int &nlocal, const int &newton_bond,\
      const KK_FLOAT &evdwl, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,\
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 protected:

  void coeff_set_tetramers_kokkos(int narg, char **arg);

  int oxdnaflag;
  enum EnabledOXDNAFlag{OXDNA=1,OXDNA3=2};

  class NeighborKokkos *neighborKK;

  t_oxdna_packed_sub<DeviceType> x;    // positions in the packed record of fix OXDNA/LRF/kk
  t_oxdna_packed<DeviceType> xn;    // the whole packed record, for row loads
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;
  typename AT::t_int_2d_lr bondlist;
  typename AT::t_tagint_1d tag;
  typename AT::t_tagint_1d id5p;
  typename AT::t_tagint_1d id3p;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int nbondlist;
  int nlocal, newton_bond, eflag, vflag;

  // stacking interaction parameters
  // all per-type-pair coefficients of a pair packed in one struct
  Kokkos::DualView<ParamsOxdnaStk2 **, Kokkos::LayoutRight, DeviceType> k_params2_st;
  typename Kokkos::DualView<ParamsOxdnaStk2 **, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params2_st;
  // all per-tetramer coefficients of a pair packed in one struct
  Kokkos::DualView<ParamsOxdnaStk4 ****, Kokkos::LayoutRight, DeviceType> k_params4_st;
  typename Kokkos::DualView<ParamsOxdnaStk4 ****, Kokkos::LayoutRight, DeviceType>::t_dev_const_randomread d_params4_st;
  // per-atom arrays for local unit vectors
  t_oxdna_packed_sub<DeviceType> d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

  void allocate() override;

  friend void pair_virial_fdotr_compute<PairOxdnaStkKokkos>(PairOxdnaStkKokkos*);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
  FixOxdnaPrimeNeighsKokkos<DeviceType> *fix_oxdna_prime_neighsKK;    // ptr to OXDNA/PRIME_NEIGHS/kk fix
  bigint last_prime_neighs_bond_nbuild;
  int tetramer_uniform;    // 1 if the tetramer coefficients do not depend on the 3'/5' context types

  // Precomputed atom a/b 3'/5' directionality and atom mapping of their 3' and 5' neighbors.
  // 0-3 : atom a, atom b, id3p[a], id5p[b] for each bond.
  typename AT::t_int_1d_4 d_prime_neighs_bond_own;
  typename AT::t_int_1d_4_randomread d_prime_neighs_bond;
};

}    // namespace LAMMPS_NS

#endif
#endif

