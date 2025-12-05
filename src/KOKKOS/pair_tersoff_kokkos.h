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
PairStyle(tersoff/kk,PairTersoffKokkos<LMPDeviceType>);
PairStyle(tersoff/kk/device,PairTersoffKokkos<LMPDeviceType>);
PairStyle(tersoff/kk/host,PairTersoffKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_TERSOFF_KOKKOS_H
#define LMP_PAIR_TERSOFF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_tersoff.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

struct ParamKokkos {
  KK_FLOAT lam1, lam2, lam3;
  KK_FLOAT c, d, h;
  KK_FLOAT gamma, powerm;
  KK_FLOAT powern, beta;
  KK_FLOAT biga, bigb, bigd, bigr;
  KK_FLOAT cut, cutsq;
  KK_FLOAT c1, c2, c3, c4;
  int ielement, jelement, kelement;
  int powermint;
  KK_FLOAT Z_i, Z_j;    // added for TersoffZBL
  KK_FLOAT ZBLcut, ZBLexpscale;
  KK_FLOAT c5, ca1, ca4;    // added for TersoffMOD
  KK_FLOAT powern_del;
  KK_FLOAT c0;    // added for TersoffMODC

  // convenient = operator
  ParamKokkos& operator=(const PairTersoff::Param& other) {
    lam1 = static_cast<KK_FLOAT>(other.lam1);
    lam2 = static_cast<KK_FLOAT>(other.lam2);
    lam3 = static_cast<KK_FLOAT>(other.lam3);
    c = static_cast<KK_FLOAT>(other.c);
    d = static_cast<KK_FLOAT>(other.d);
    h = static_cast<KK_FLOAT>(other.h);
    gamma = static_cast<KK_FLOAT>(other.gamma);
    powerm = static_cast<KK_FLOAT>(other.powerm);
    powern = static_cast<KK_FLOAT>(other.powern);
    beta = static_cast<KK_FLOAT>(other.beta);
    biga = static_cast<KK_FLOAT>(other.biga);
    bigb = static_cast<KK_FLOAT>(other.bigb);
    bigd = static_cast<KK_FLOAT>(other.bigd);
    bigr = static_cast<KK_FLOAT>(other.bigr);
    cut = static_cast<KK_FLOAT>(other.cut);
    cutsq = static_cast<KK_FLOAT>(other.cutsq);
    c1 = static_cast<KK_FLOAT>(other.c1);
    c2 = static_cast<KK_FLOAT>(other.c2);
    c3 = static_cast<KK_FLOAT>(other.c3);
    c4 = static_cast<KK_FLOAT>(other.c4);
    ielement = other.ielement;
    jelement = other.jelement;
    kelement = other.kelement;
    powermint = other.powermint;
    Z_i = static_cast<KK_FLOAT>(other.Z_i);
    Z_j = static_cast<KK_FLOAT>(other.Z_j);
    ZBLcut = static_cast<KK_FLOAT>(other.ZBLcut);
    ZBLexpscale = static_cast<KK_FLOAT>(other.ZBLexpscale);
    c5 = static_cast<KK_FLOAT>(other.c5);
    ca1 = static_cast<KK_FLOAT>(other.ca1);
    ca4 = static_cast<KK_FLOAT>(other.ca4);
    powern_del = static_cast<KK_FLOAT>(other.powern_del);
    c0 = static_cast<KK_FLOAT>(other.c0);
    return *this;
  }

};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairTersoffCompute{};

struct TagPairTersoffComputeShortNeigh{};

template<class DeviceType>
class PairTersoffKokkos : public PairTersoff {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  // Static blocking size for PairTersoffCompute, EVFLAG == 0
  static constexpr int block_size_compute_tersoff_force = 128;
  // EVFLAG == 1, intentionally different due to how Kokkos implements
  // reductions vs simple parallel_for
  static constexpr int block_size_compute_tersoff_energy = 256;

  PairTersoffKokkos(class LAMMPS *);
  ~PairTersoffKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  // RangePolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&) const;

  // MDRangePolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const int&, const int&) const;

  // TeamPolicy versions
  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagPairTersoffCompute<NEIGHFLAG,EVFLAG> >::member_type&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffCompute<NEIGHFLAG,EVFLAG>, const typename Kokkos::TeamPolicy<DeviceType, TagPairTersoffCompute<NEIGHFLAG,EVFLAG> >::member_type&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairTersoffComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void tersoff_compute(const int&, EV_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_fc_k(const ParamKokkos& param, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dfc(const ParamKokkos& param, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  void ters_fc_k_and_ters_dfc(const ParamKokkos& param, const KK_FLOAT &r, KK_FLOAT &fc, KK_FLOAT &dfc) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_fa_k(const ParamKokkos& param, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dfa(const ParamKokkos& param, const KK_FLOAT &r) const;

  KOKKOS_INLINE_FUNCTION
  void ters_fa_k_and_ters_dfa(const ParamKokkos& param, const KK_FLOAT &r, KK_FLOAT &fa, KK_FLOAT &dfa) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_bij_k(const ParamKokkos& param, const KK_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dbij(const ParamKokkos& param, const KK_FLOAT &bo) const;

  KOKKOS_INLINE_FUNCTION
  void ters_bij_k_and_ters_dbij(const ParamKokkos& param, const KK_FLOAT &bo, KK_FLOAT &bij, KK_FLOAT &prefactor) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT bondorder(const ParamKokkos& param,
              const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
              const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_gijk(const ParamKokkos& param, const KK_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT ters_dgijk(const ParamKokkos& param, const KK_FLOAT &cos) const;

  KOKKOS_INLINE_FUNCTION
  void ters_gijk_and_ters_dgijk(const ParamKokkos& param, const KK_FLOAT &cos, KK_FLOAT& gijk, KK_FLOAT& dgijk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthb(const ParamKokkos& param, const KK_FLOAT &prefactor,
              const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
              const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
              KK_ACC_FLOAT *fi, KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbj(const ParamKokkos& param, const KK_FLOAT &prefactor,
              const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
              const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
              KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  void ters_dthbk(const ParamKokkos& param, const KK_FLOAT &prefactor,
              const KK_FLOAT &rij, const KK_FLOAT &dx1, const KK_FLOAT &dy1, const KK_FLOAT &dz1,
              const KK_FLOAT &rik, const KK_FLOAT &dx2, const KK_FLOAT &dy2, const KK_FLOAT &dz2,
              KK_ACC_FLOAT *fk) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT vec3_dot(const KK_FLOAT x[3], const KK_FLOAT y[3]) const {
    KK_FLOAT dot = 0;
    for (int i = 0; i < 3; i++)
      dot += x[i]*y[i];
    return dot;
  }

  KOKKOS_INLINE_FUNCTION
  void vec3_add(const KK_FLOAT x[3], const KK_FLOAT y[3], KK_FLOAT * const z) const {
    for (int i = 0; i < 3; i++)
      z[i] = x[i]+y[i];
  }

  template<typename k_type, typename x_type, typename y_type>
  KOKKOS_INLINE_FUNCTION
  void vec3_scale(const k_type k, const x_type x[3], y_type y[3]) const {
    for (int i = 0; i < 3; i++)
      y[i] = static_cast<y_type>(static_cast<x_type>(k)*x[i]);
  }

  template<typename kx_type, typename yz_type>
  KOKKOS_INLINE_FUNCTION
  void vec3_scaleadd(const kx_type k, const kx_type x[3], const yz_type y[3], yz_type z[3]) const {
    for (int i = 0; i < 3; i++)
      z[i] = static_cast<yz_type>(k*x[i])+y[i];
  }

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void v_tally3(EV_FLOAT &ev, const int &i, const int &j, const int &k,
                KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drij, KK_FLOAT *drik) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally3_atom(EV_FLOAT &ev, const int &i, const int &j, const int &k,
                KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drjk) const;

  void setup_params() override;

 protected:
  typename AT::t_int_3d_randomread d_elem3param;
  typename AT::t_int_1d_randomread d_map;

  typedef Kokkos::DualView<ParamKokkos*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  KK_FLOAT cutmax_sq;

  int inum;
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int need_dup;

  using KKDeviceType = typename KKDevice<DeviceType>::value;

  template<typename DataType, typename Layout>
  using DupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterDuplicated>;

  template<typename DataType, typename Layout>
  using NonDupScatterView = KKScatterView<DataType, Layout, KKDeviceType, KKScatterSum, KKScatterNonDuplicated>;

  DupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> dup_f;
  DupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> dup_eatom;
  DupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> dup_vatom;

  NonDupScatterView<KK_ACC_FLOAT*[3], typename DAT::t_kkacc_1d_3::array_layout> ndup_f;
  NonDupScatterView<KK_ACC_FLOAT*, typename DAT::t_kkacc_1d::array_layout> ndup_eatom;
  NonDupScatterView<KK_ACC_FLOAT*[6], typename DAT::t_kkacc_1d_6::array_layout> ndup_vatom;

  typedef Kokkos::DualView<KK_FLOAT**[7],Kokkos::LayoutRight,DeviceType> tdual_kkfloat_2d_n7;
  typedef typename tdual_kkfloat_2d_n7::t_dev_const_randomread t_kkfloat_2d_n7_randomread;
  typedef typename tdual_kkfloat_2d_n7::t_host t_hostkkfloat_2d_n7;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<DeviceType> k_list;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  typename AT::t_int_2d_dl d_neighbors_short;
  typename AT::t_int_1d d_numneigh_short;

  friend void pair_virial_fdotr_compute<PairTersoffKokkos>(PairTersoffKokkos*);
};

}

#endif
#endif

