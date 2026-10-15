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

/* ----------------------------------------------------------------------
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(dpd/coul/slater/long/kk,PairDPDCoulSlaterLongKokkos<LMPDeviceType>);
PairStyle(dpd/coul/slater/long/kk/device,PairDPDCoulSlaterLongKokkos<LMPDeviceType>);
PairStyle(dpd/coul/slater/long/kk/host,PairDPDCoulSlaterLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_DPD_COUL_SLATER_LONG_KOKKOS_H
#define LMP_PAIR_DPD_COUL_SLATER_LONG_KOKKOS_H

#include "pair_dpd_coul_slater_long.h"
#include "pair_kokkos.h"
#include "kokkos_type.h"

// a build configured with -D KOKKOS_DEBUG_RNG=on draws the random numbers from
// the host RanMars generator, so that the results can be compared bit for bit
// with the plain style.  an explicit choice of generator still wins.

#if defined(LMP_KOKKOS_DEBUG_RNG) && !defined(DPD_USE_Random_XorShift64) && !defined(Random_XorShift1024)
#define DPD_USE_RAN_MARS
#endif

#if !defined(DPD_USE_RAN_MARS) && !defined(DPD_USE_Random_XorShift64) && !defined(Random_XorShift1024)
#define DPD_USE_Random_XorShift64
#endif

#ifdef DPD_USE_RAN_MARS
#include "rand_pool_wrap_kokkos.h"
#else
#include "Kokkos_Random.hpp"
#endif

namespace LAMMPS_NS {

template<class DeviceType>
class PairDPDCoulSlaterLongKokkos : public PairDPDCoulSlaterLong {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  PairDPDCoulSlaterLongKokkos(class LAMMPS*);
  ~PairDPDCoulSlaterLongKokkos() override;

  void allocate() override;

  void init_style() override;
  double init_one(int i, int j) override;
  void compute(int, int) override;

  struct params_dpd {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_dpd() {cut_dpd=cut_dpdsq=cut_slatersq=a0=gamma=sigma=0;}
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_dpd(int /*i*/) {cut_dpd=cut_dpdsq=cut_slatersq=a0=gamma=sigma=0;}
    KK_FLOAT cut_dpd,cut_dpdsq,cut_slatersq,a0,gamma,sigma;
  };

  template<int NEIGHFLAG, int EVFLAG>
  struct TagDPDCoulSlaterLongKokkos{};

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator () (TagDPDCoulSlaterLongKokkos<NEIGHFLAG,EVFLAG>, const int &i) const;

  template<int NEIGHFLAG, int EVFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator () (TagDPDCoulSlaterLongKokkos<NEIGHFLAG,EVFLAG>, const int &i, EV_FLOAT&) const;

  template<int NEIGHFLAG>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
                const KK_FLOAT &epair, const KK_FLOAT &fpair,
                const KK_FLOAT &delx, const KK_FLOAT &dely, const KK_FLOAT &delz) const;

 private:
  KK_FLOAT special_lj[4], special_rf[4], special_coul[4];
  int eflag,vflag;
  int neighflag,nlocal;
  KK_FLOAT dtinvsqrt;
  KK_FLOAT qqrd2e, g_ewald_kk, lamda_kk;

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

#ifdef DPD_USE_RAN_MARS
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#elif defined(DPD_USE_Random_XorShift64)
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#elif defined(DPD_USE_Random_XorShift1024)
  Kokkos::Random_XorShift1024_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift1024_Pool<DeviceType>::generator_type rand_type;
#endif
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_randomread v;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_kkfloat_1d_randomread q;
  typename AT::t_int_1d_randomread type;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  Kokkos::DualView<params_dpd**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_dpd**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;
  friend void pair_virial_fdotr_compute<PairDPDCoulSlaterLongKokkos>(PairDPDCoulSlaterLongKokkos*);

};
}
#endif
#endif
