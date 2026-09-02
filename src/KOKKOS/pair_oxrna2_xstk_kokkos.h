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
PairStyle(oxrna2/xstk/kk,PairOxrna2XstkKokkos<LMPDeviceType>);
PairStyle(oxrna2/xstk/kk/device,PairOxrna2XstkKokkos<LMPDeviceType>);
PairStyle(oxrna2/xstk/kk/host,PairOxrna2XstkKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_OXRNA2_XSTK_KOKKOS_H
#define LMP_PAIR_OXRNA2_XSTK_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_oxrna2_xstk.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxrna2XstkCompute {};

template<class DeviceType>
class PairOxrna2XstkKokkos : public PairOxrna2Xstk, public KokkosBase {
 public:
  enum {EnabledNeighFlags = FULL | HALFTHREAD | HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairOxrna2XstkKokkos(class LAMMPS *);
  ~PairOxrna2XstkKokkos() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxrna2XstkCompute<NEIGHFLAG, NEWTON_PAIR, EVFLAG>,
                  const int &, EV_FLOAT &) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxrna2XstkCompute<NEIGHFLAG, NEWTON_PAIR, EVFLAG>,
                  const int &) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j, const KK_FLOAT &epair,
                    const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy,
                    const KK_ACC_FLOAT &fz, const KK_FLOAT &delx,
                    const KK_FLOAT &dely, const KK_FLOAT &delz) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int &j) const;

 protected:
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

  int neighflag;
  int nlocal, eflag, vflag;
  int anum;

  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;

  // cross-stacking interaction parameters
  typename AT::tdual_kkfloat_2d k_k_xst, k_cut_xst_0, k_cut_xst_c;
  typename AT::tdual_kkfloat_2d k_cut_xst_lo, k_cut_xst_hi;
  typename AT::tdual_kkfloat_2d k_cut_xst_lc, k_cut_xst_hc, k_b_xst_lo, k_b_xst_hi;
  typename AT::tdual_kkfloat_2d k_cutsq_xst_hc;
  typename AT::tdual_kkfloat_2d k_a_xst1, k_theta_xst1_0, k_dtheta_xst1_ast;
  typename AT::tdual_kkfloat_2d k_b_xst1, k_dtheta_xst1_c;
  typename AT::tdual_kkfloat_2d k_a_xst2, k_theta_xst2_0, k_dtheta_xst2_ast;
  typename AT::tdual_kkfloat_2d k_b_xst2, k_dtheta_xst2_c;
  typename AT::tdual_kkfloat_2d k_a_xst3, k_theta_xst3_0, k_dtheta_xst3_ast;
  typename AT::tdual_kkfloat_2d k_b_xst3, k_dtheta_xst3_c;
  typename AT::tdual_kkfloat_2d k_a_xst7, k_theta_xst7_0, k_dtheta_xst7_ast;
  typename AT::tdual_kkfloat_2d k_b_xst7, k_dtheta_xst7_c;
  typename AT::tdual_kkfloat_2d k_a_xst8, k_theta_xst8_0, k_dtheta_xst8_ast;
  typename AT::tdual_kkfloat_2d k_b_xst8, k_dtheta_xst8_c;
  typename AT::t_kkfloat_2d_randomread d_k_xst, d_cut_xst_0, d_cut_xst_c;
  typename AT::t_kkfloat_2d_randomread d_cut_xst_lo, d_cut_xst_hi;
  typename AT::t_kkfloat_2d_randomread d_cut_xst_lc, d_cut_xst_hc, d_b_xst_lo, d_b_xst_hi;
  typename AT::t_kkfloat_2d_randomread d_cutsq_xst_hc;
  typename AT::t_kkfloat_2d_randomread d_a_xst1, d_theta_xst1_0, d_dtheta_xst1_ast;
  typename AT::t_kkfloat_2d_randomread d_b_xst1, d_dtheta_xst1_c;
  typename AT::t_kkfloat_2d_randomread d_a_xst2, d_theta_xst2_0, d_dtheta_xst2_ast;
  typename AT::t_kkfloat_2d_randomread d_b_xst2, d_dtheta_xst2_c;
  typename AT::t_kkfloat_2d_randomread d_a_xst3, d_theta_xst3_0, d_dtheta_xst3_ast;
  typename AT::t_kkfloat_2d_randomread d_b_xst3, d_dtheta_xst3_c;
  typename AT::t_kkfloat_2d_randomread d_a_xst7, d_theta_xst7_0, d_dtheta_xst7_ast;
  typename AT::t_kkfloat_2d_randomread d_b_xst7, d_dtheta_xst7_c;
  typename AT::t_kkfloat_2d_randomread d_a_xst8, d_theta_xst8_0, d_dtheta_xst8_ast;
  typename AT::t_kkfloat_2d_randomread d_b_xst8, d_dtheta_xst8_c;

  // per-atom arrays for local unit vectors
  typename AT::t_kkfloat_1d_3_randomread d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, \
  KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*, typename AT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename AT::t_kkacc_1d_6::array_layout> dup_vatom;
  NonDupScatterView<KK_ACC_FLOAT*, typename AT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename AT::t_kkacc_1d_6::array_layout> ndup_vatom;

  void allocate() override;

  friend void pair_virial_fdotr_compute<PairOxrna2XstkKokkos>(PairOxrna2XstkKokkos *);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;  // ptr to OXDNA/LRF/kk fix
};

}    // namespace LAMMPS_NS

#endif
#endif
