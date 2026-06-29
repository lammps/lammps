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
PairStyle(lj/switch3/coulgauss/long/kk,PairLJSwitch3CoulGaussLongKokkos<LMPDeviceType>);
PairStyle(lj/switch3/coulgauss/long/kk/device,PairLJSwitch3CoulGaussLongKokkos<LMPDeviceType>);
PairStyle(lj/switch3/coulgauss/long/kk/host,PairLJSwitch3CoulGaussLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_SWITCH3_COULGAUSS_LONG_KOKKOS_H
#define LMP_PAIR_LJ_SWITCH3_COULGAUSS_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_switch3_coulgauss_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJSwitch3CoulGaussLongKokkos : public PairLJSwitch3CoulGaussLong {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJSwitch3CoulGaussLongKokkos(class LAMMPS *);
  ~PairLJSwitch3CoulGaussLongKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_lj_sw3_cg {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_sw3_cg() {cutsq=0;lj2=0;lj3=0;lj4=0;offset=0;
                        cut_lj=0;cut_ljsq=0;cut_coulsq=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_sw3_cg(int /*i*/) {cutsq=0;lj2=0;lj3=0;lj4=0;offset=0;
                                  cut_lj=0;cut_ljsq=0;cut_coulsq=0;};
    KK_FLOAT cutsq,lj2,lj3,lj4,offset,cut_lj,cut_ljsq,cut_coulsq;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int& j,
                         const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int& j,
                         const int& itype, const int& jtype,
                         const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int& j,
                          const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int& j,
                          const int& itype, const int& jtype,
                          const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_sw3_cg**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_sw3_cg**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj_sw3_cg m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  DAT::ttransform_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT g_ewald_kk;
  KK_FLOAT truncw_kk;
  KK_FLOAT truncwi_kk;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJSwitch3CoulGaussLongKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJSwitch3CoulGaussLongKokkos,FULL,0>(PairLJSwitch3CoulGaussLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSwitch3CoulGaussLongKokkos,FULL,1>(PairLJSwitch3CoulGaussLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSwitch3CoulGaussLongKokkos,HALF>(PairLJSwitch3CoulGaussLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSwitch3CoulGaussLongKokkos,HALFTHREAD>(PairLJSwitch3CoulGaussLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJSwitch3CoulGaussLongKokkos,void>(PairLJSwitch3CoulGaussLongKokkos*,
                                                                       NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJSwitch3CoulGaussLongKokkos>(PairLJSwitch3CoulGaussLongKokkos*);
};

}

#endif
#endif
