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
PairStyle(edip/kk,PairEDIPKokkos<LMPDeviceType>);
PairStyle(edip/kk/device,PairEDIPKokkos<LMPDeviceType>);
PairStyle(edip/kk/host,PairEDIPKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_EDIP_KOKKOS_H
#define LMP_PAIR_EDIP_KOKKOS_H

#include "pair_edip.h"
#include "pair_kokkos.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairEDIPCompute{};

struct TagPairEDIPPreLoop{};

namespace LAMMPS_NS {

template<class DeviceType>
class PairEDIPKokkos : public PairEDIP {
 public:
  enum {EnabledNeighFlags=HALF|HALFTHREAD};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairEDIPKokkos(class LAMMPS *);
  ~PairEDIPKokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEDIPCompute<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEDIPCompute<NEIGHFLAG,EVFLAG>, const int&) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairEDIPPreLoop, const int&) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
            const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                       KK_ACC_FLOAT *fj, KK_ACC_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

 protected:
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT interpolate(const typename AT::t_kkfloat_1d_randomread &, const int &, const KK_FLOAT &) const;

  // copies of the interpolation grids of the base class

  DAT::tdual_kkfloat_1d k_cutoffFunction, k_cutoffFunctionDerived, k_pow2B, k_exp2B, k_exp3B;
  DAT::tdual_kkfloat_1d k_qFunctionGrid, k_expMinusBetaZeta_iZeta_iGrid;
  DAT::tdual_kkfloat_1d k_tauFunctionGrid, k_tauFunctionDerivedGrid;
  typename AT::t_kkfloat_1d_randomread d_cutoffFunction, d_cutoffFunctionDerived, d_pow2B, d_exp2B, d_exp3B;
  typename AT::t_kkfloat_1d_randomread d_qFunctionGrid, d_expMinusBetaZeta_iZeta_iGrid;
  typename AT::t_kkfloat_1d_randomread d_tauFunctionGrid, d_tauFunctionDerivedGrid;

  KK_FLOAT cutsq_kk;

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;

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

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  int inum;

  // neighbors within the cutoff with the per-pair quantities of the three-body loop,
  // and the coordination Z of each atom

  typename AT::t_int_2d_dl d_neighbors_short;
  typename AT::t_int_1d d_numneigh_short;
  typename AT::t_kkfloat_2d_dl d_preInvR, d_preExp3B, d_preExp3BDerived;
  typename AT::t_kkfloat_1d d_zeta;

  DAT::tdual_int_scalar k_overflow;
  typename AT::t_int_scalar d_overflow;

  friend void pair_virial_fdotr_compute<PairEDIPKokkos>(PairEDIPKokkos*);
};

}

#endif
#endif
