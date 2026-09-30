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
PairStyle(lj/smooth/kk,PairLJSmoothKokkos<LMPDeviceType>);
PairStyle(lj/smooth/kk/device,PairLJSmoothKokkos<LMPDeviceType>);
PairStyle(lj/smooth/kk/host,PairLJSmoothKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_SMOOTH_KOKKOS_H
#define LMP_PAIR_LJ_SMOOTH_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_smooth.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJSmoothKokkos : public PairLJSmooth {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJSmoothKokkos(class LAMMPS *);
  ~PairLJSmoothKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj_smooth{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_smooth() {cutsq=0;cut_inner_sq=0;lj1=0;lj2=0;lj3=0;lj4=0;
                        ljsw0=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;
                        cut_inner=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_smooth(int /*i*/) {cutsq=0;cut_inner_sq=0;lj1=0;lj2=0;lj3=0;lj4=0;
                        ljsw0=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;
                        cut_inner=0;offset=0;};
    KK_FLOAT cutsq,cut_inner_sq,lj1,lj2,lj3,lj4;
    KK_FLOAT ljsw0,ljsw1,ljsw2,ljsw3,ljsw4;
    KK_FLOAT cut_inner,offset;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  Kokkos::DualView<params_lj_smooth**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_smooth**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj_smooth m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;
  KK_FLOAT special_lj[4];

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJSmoothKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJSmoothKokkos,FULL,0>(PairLJSmoothKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSmoothKokkos,FULL,1>(PairLJSmoothKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSmoothKokkos,HALF>(PairLJSmoothKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSmoothKokkos,HALFTHREAD>(PairLJSmoothKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJSmoothKokkos>(PairLJSmoothKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJSmoothKokkos>(PairLJSmoothKokkos*);
};

}

#endif
#endif
