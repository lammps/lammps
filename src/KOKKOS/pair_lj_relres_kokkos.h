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
PairStyle(lj/relres/kk,PairLJRelResKokkos<LMPDeviceType>);
PairStyle(lj/relres/kk/device,PairLJRelResKokkos<LMPDeviceType>);
PairStyle(lj/relres/kk/host,PairLJRelResKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_RELRES_KOKKOS_H
#define LMP_PAIR_LJ_RELRES_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_relres.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJRelResKokkos : public PairLJRelRes {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJRelResKokkos(class LAMMPS *);
  ~PairLJRelResKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj() {cutsq=0;cutf_inner_sq=0;cutf_inner=0;cutfsq=0;cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;ljf1=0;ljf2=0;ljf3=0;ljf4=0;ljsw0=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;ljswf0=0;ljswf1=0;ljswf2=0;ljswf3=0;ljswf4=0;offset=0;offsetsp=0;offsetsm=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj(int /*i*/) {cutsq=0;cutf_inner_sq=0;cutf_inner=0;cutfsq=0;cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;ljf1=0;ljf2=0;ljf3=0;ljf4=0;ljsw0=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;ljswf0=0;ljswf1=0;ljswf2=0;ljswf3=0;ljswf4=0;offset=0;offsetsp=0;offsetsm=0;};
    KK_FLOAT cutsq,cutf_inner_sq,cutf_inner,cutfsq,cut_inner_sq,cut_inner,lj1,lj2,lj3,lj4,ljf1,ljf2,ljf3,ljf4,ljsw0,ljsw1,ljsw2,ljsw3,ljsw4,ljswf0,ljswf1,ljswf2,ljswf3,ljswf4,offset,offsetsp,offsetsm;
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

  Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];  // hardwired to space for 12 atom types
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
  friend struct PairComputeFunctor<PairLJRelResKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJRelResKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJRelResKokkos,FULL,0>(PairLJRelResKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJRelResKokkos,FULL,1>(PairLJRelResKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJRelResKokkos,HALF>(PairLJRelResKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJRelResKokkos,HALFTHREAD>(PairLJRelResKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJRelResKokkos>(PairLJRelResKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJRelResKokkos>(PairLJRelResKokkos*);
};

}

#endif
#endif

