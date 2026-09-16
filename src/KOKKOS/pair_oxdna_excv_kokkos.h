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
PairStyle(oxdna/excv/kk,PairOxdnaExcvKokkos<LMPDeviceType>);
PairStyle(oxdna/excv/kk/device,PairOxdnaExcvKokkos<LMPDeviceType>);
PairStyle(oxdna/excv/kk/host,PairOxdnaExcvKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA_EXCV_KOKKOS_H
#define LMP_PAIR_OXDNA_EXCV_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_oxdna_excv.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<class DeviceType>
class FixOxdnaNpairKokkos;  // forward declaration

template<class DeviceType>
class FixOxdnaPrimeNeighsKokkos;  // forward declaration

template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaExcvCompute{};

template<int NEIGHFLAG, int NEWTON_PAIR>
struct ev_tally_xyz{};

template<class DeviceType>
class PairOxdnaExcvKokkos : public PairOxdnaExcv, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdnaExcvKokkos(class LAMMPS *);
  ~PairOxdnaExcvKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void coeff(int, char **) override;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int OXDNAFLAG, int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaExcvCompute<OXDNAFLAG,NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:

  void coeff_set_tetramers_kokkos(int narg, char **arg);

  int oxdnaflag;
  enum EnabledOXDNAFlag{OXDNA=1,OXDNA2=2,OXRNA2=4,OXDNA3=8};

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkacc_1d_3 torque;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;
  double special_lj[4];

  typename AT::tdual_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int neighflag;
  int nlocal, eflag, vflag;
  int anum;
  bigint last_prime_neighs_pair_lastcall;

  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;
  typename AT::t_int_3d_randomread d_prime_neighs_pair;

  // s=sugar-phosphate backbone site, b=base site, st=stacking site
  // excluded volume interaction parameters
  typename AT::tdual_kkfloat_2d k_epsilon_bkbk, k_sigma_bkbk, k_cut_bkbk_ast, k_cutsq_bkbk_ast;
  typename AT::tdual_kkfloat_2d k_lj1_bkbk, k_lj2_bkbk, k_b_bkbk, k_cut_bkbk_c, k_cutsq_bkbk_c;
  typename AT::tdual_kkfloat_2d k_epsilon_bkbs, k_sigma_bkbs, k_cut_bkbs_ast, k_cutsq_bkbs_ast;
  typename AT::tdual_kkfloat_2d k_lj1_bkbs, k_lj2_bkbs, k_b_bkbs, k_cut_bkbs_c, k_cutsq_bkbs_c;
  typename AT::tdual_kkfloat_2d k_epsilon_bsbs, k_sigma_bsbs, k_cut_bsbs_ast, k_cutsq_bsbs_ast;
  typename AT::tdual_kkfloat_2d k_lj1_bsbs, k_lj2_bsbs, k_b_bsbs, k_cut_bsbs_c, k_cutsq_bsbs_c;
  typename AT::t_kkfloat_2d_randomread d_epsilon_bkbk, d_sigma_bkbk, d_cut_bkbk_ast, d_cutsq_bkbk_ast;
  typename AT::t_kkfloat_2d_randomread d_lj1_bkbk, d_lj2_bkbk, d_b_bkbk, d_cut_bkbk_c, d_cutsq_bkbk_c;
  typename AT::t_kkfloat_2d_randomread d_epsilon_bkbs, d_sigma_bkbs, d_cut_bkbs_ast, d_cutsq_bkbs_ast;
  typename AT::t_kkfloat_2d_randomread d_lj1_bkbs, d_lj2_bkbs, d_b_bkbs, d_cut_bkbs_c, d_cutsq_bkbs_c;
  typename AT::t_kkfloat_2d_randomread d_epsilon_bsbs, d_sigma_bsbs, d_cut_bsbs_ast, d_cutsq_bsbs_ast;
  typename AT::t_kkfloat_2d_randomread d_lj1_bsbs, d_lj2_bsbs, d_b_bsbs, d_cut_bsbs_c, d_cutsq_bsbs_c;
  // tetramer-dependent coefficients
  typename AT::tdual_kkfloat_4d k_sigma4_bsbs, k_cut4_bsbs_ast, k_cut4sq_bsbs_ast;
  typename AT::tdual_kkfloat_4d k_lj14_bsbs, k_lj24_bsbs, k_b4_bsbs, k_cut4_bsbs_c, k_cut4sq_bsbs_c;
  typename AT::t_kkfloat_4d_randomread d_sigma4_bsbs, d_cut4_bsbs_ast, d_cut4sq_bsbs_ast;
  typename AT::t_kkfloat_4d_randomread d_lj14_bsbs, d_lj24_bsbs, d_b4_bsbs, d_cut4_bsbs_c, d_cut4sq_bsbs_c;
  // per-atom arrays for local unit vectors
  DAT::tdual_kkfloat_1d_3 k_nx_xtrct, k_ny_xtrct, k_nz_xtrct;
  typename AT::t_kkfloat_1d_3_randomread d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread tag;
  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread id5p;
  typename ArrayTypes<DeviceType>::t_tagint_1d_randomread id3p;

  int map_style;
  DAT::tdual_int_1d k_map_array;
  dual_hash_type k_map_hash;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> dup_torque;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*[3], typename AT::t_kkacc_1d_3::array_layout> ndup_torque;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  void allocate() override;

  friend void pair_virial_fdotr_compute<PairOxdnaExcvKokkos>(PairOxdnaExcvKokkos*);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
  FixOxdnaNpairKokkos<DeviceType> *fix_oxdna_npairKK;    // ptr to OXDNA/NPAIR/kk fix
  FixOxdnaPrimeNeighsKokkos<DeviceType> *fix_oxdna_prime_neighsKK;    // ptr to OXDNA/PRIME_NEIGHS/kk fix
};

}

#endif
#endif

