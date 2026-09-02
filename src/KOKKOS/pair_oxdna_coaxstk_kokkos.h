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
PairStyle(oxdna/coaxstk/kk,PairOxdnaCoaxstkKokkos<LMPDeviceType>);
PairStyle(oxdna/coaxstk/kk/device,PairOxdnaCoaxstkKokkos<LMPDeviceType>);
PairStyle(oxdna/coaxstk/kk/host,PairOxdnaCoaxstkKokkos<LMPHostType>);
PairStyle(oxrna2/coaxstk/kk,PairOxdnaCoaxstkKokkos<LMPDeviceType>);
PairStyle(oxrna2/coaxstk/kk/device,PairOxdnaCoaxstkKokkos<LMPDeviceType>);
PairStyle(oxrna2/coaxstk/kk/host,PairOxdnaCoaxstkKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_OXDNA_COAXSTK_KOKKOS_H
#define LMP_PAIR_OXDNA_COAXSTK_KOKKOS_H

#include "kokkos_base.h"
#include "pair_kokkos.h"
#include "pair_oxdna_coaxstk.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixOxdnaLRFKokkos;  // forward declaration

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairOxdnaCoaxstkCompute{};

template<class DeviceType>
class PairOxdnaCoaxstkKokkos : public PairOxdnaCoaxstk, public KokkosBase {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairOxdnaCoaxstkKokkos(class LAMMPS *);
  ~PairOxdnaCoaxstkKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaCoaxstkCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairOxdnaCoaxstkCompute<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally_xyz(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_ACC_FLOAT &fx, const KK_ACC_FLOAT &fy, const KK_ACC_FLOAT &fz,
      const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

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

  typename AT::tdual_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int neighflag;
  int nlocal, eflag, vflag;
  int anum;

  typename AT::t_neighbors_2d_randomread d_neighbors;
  typename AT::t_int_1d_randomread d_alist;
  typename AT::t_int_1d_randomread d_numneigh;

  // coaxial stacking interaction parameters
  typename AT::tdual_kkfloat_2d k_k_cxst, k_cut_cxst_0, k_cut_cxst_c;
  typename AT::tdual_kkfloat_2d k_cut_cxst_lo, k_cut_cxst_hi;
  typename AT::tdual_kkfloat_2d k_cut_cxst_lc, k_cut_cxst_hc, k_b_cxst_lo, k_b_cxst_hi;
  typename AT::tdual_kkfloat_2d k_cutsq_cxst_hc;
  typename AT::tdual_kkfloat_2d k_a_cxst1, k_theta_cxst1_0, k_dtheta_cxst1_ast;
  typename AT::tdual_kkfloat_2d k_b_cxst1, k_dtheta_cxst1_c;
  typename AT::tdual_kkfloat_2d k_a_cxst4, k_theta_cxst4_0, k_dtheta_cxst4_ast;
  typename AT::tdual_kkfloat_2d k_b_cxst4, k_dtheta_cxst4_c;
  typename AT::tdual_kkfloat_2d k_a_cxst5, k_theta_cxst5_0, k_dtheta_cxst5_ast;
  typename AT::tdual_kkfloat_2d k_b_cxst5, k_dtheta_cxst5_c;
  typename AT::tdual_kkfloat_2d k_a_cxst6, k_theta_cxst6_0, k_dtheta_cxst6_ast;
  typename AT::tdual_kkfloat_2d k_b_cxst6, k_dtheta_cxst6_c;
  typename AT::tdual_kkfloat_2d k_a_cxst3p, k_cosphi_cxst3p_ast, k_b_cxst3p, k_cosphi_cxst3p_c;
  typename AT::tdual_kkfloat_2d k_a_cxst4p, k_cosphi_cxst4p_ast, k_b_cxst4p, k_cosphi_cxst4p_c;
  typename AT::t_kkfloat_2d_randomread d_k_cxst, d_cut_cxst_0, d_cut_cxst_c;
  typename AT::t_kkfloat_2d_randomread d_cut_cxst_lo, d_cut_cxst_hi;
  typename AT::t_kkfloat_2d_randomread d_cut_cxst_lc, d_cut_cxst_hc, d_b_cxst_lo, d_b_cxst_hi;
  typename AT::t_kkfloat_2d_randomread d_cutsq_cxst_hc;
  typename AT::t_kkfloat_2d_randomread d_a_cxst1, d_theta_cxst1_0, d_dtheta_cxst1_ast;
  typename AT::t_kkfloat_2d_randomread d_b_cxst1, d_dtheta_cxst1_c;
  typename AT::t_kkfloat_2d_randomread d_a_cxst4, d_theta_cxst4_0, d_dtheta_cxst4_ast;
  typename AT::t_kkfloat_2d_randomread d_b_cxst4, d_dtheta_cxst4_c;
  typename AT::t_kkfloat_2d_randomread d_a_cxst5, d_theta_cxst5_0, d_dtheta_cxst5_ast;
  typename AT::t_kkfloat_2d_randomread d_b_cxst5, d_dtheta_cxst5_c;
  typename AT::t_kkfloat_2d_randomread d_a_cxst6, d_theta_cxst6_0, d_dtheta_cxst6_ast;
  typename AT::t_kkfloat_2d_randomread d_b_cxst6, d_dtheta_cxst6_c;
  typename AT::t_kkfloat_2d_randomread d_a_cxst3p, d_cosphi_cxst3p_ast, d_b_cxst3p, d_cosphi_cxst3p_c;
  typename AT::t_kkfloat_2d_randomread d_a_cxst4p, d_cosphi_cxst4p_ast, d_b_cxst4p, d_cosphi_cxst4p_c;
  // per-atom arrays for local unit vectors
  DAT::tdual_kkfloat_1d_3 k_nx_xtrct, k_ny_xtrct, k_nz_xtrct;
  typename AT::t_kkfloat_1d_3_randomread d_nx_xtrct, d_ny_xtrct, d_nz_xtrct;

  int first;
  typename AT::t_int_1d d_sendlist;
  typename AT::t_double_1d_um v_buf;

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

  friend void pair_virial_fdotr_compute<PairOxdnaCoaxstkKokkos>(PairOxdnaCoaxstkKokkos*);

  FixOxdnaLRFKokkos<DeviceType> *fix_oxdna_lrfKK;    // ptr to OXDNA/LRF/kk fix
};

}

#endif
#endif

