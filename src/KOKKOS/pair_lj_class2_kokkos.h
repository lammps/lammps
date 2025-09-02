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
PairStyle(lj/class2/kk,PairLJClass2Kokkos<LMPDeviceType>);
PairStyle(lj/class2/kk/device,PairLJClass2Kokkos<LMPDeviceType>);
PairStyle(lj/class2/kk/host,PairLJClass2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CLASS2_KOKKOS_H
#define LMP_PAIR_LJ_CLASS2_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_class2.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJClass2Kokkos : public PairLJClass2 {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJClass2Kokkos(class LAMMPS *);
  ~PairLJClass2Kokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_lj{
    KOKKOS_INLINE_FUNCTION
    params_lj() {cutsq=0,lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_lj(int /*i*/) {cutsq=0,lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
    KK_FLOAT cutsq,lj1,lj2,lj3,lj4,offset;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/, const int& /*itype*/,
                        const int& /*jtype*/, const KK_FLOAT& /*factor_coul*/, const KK_FLOAT& /*qtmp*/) const { return 0; }

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
  friend struct PairComputeFunctor<PairLJClass2Kokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJClass2Kokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2Kokkos,FULL,0>(PairLJClass2Kokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2Kokkos,FULL,1>(PairLJClass2Kokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2Kokkos,HALF>(PairLJClass2Kokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2Kokkos,HALFTHREAD>(PairLJClass2Kokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJClass2Kokkos>(PairLJClass2Kokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJClass2Kokkos>(PairLJClass2Kokkos*);
};

}

#endif
#endif

