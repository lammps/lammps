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
PairStyle(yukawa/colloid/kk,PairYukawaColloidKokkos<LMPDeviceType>);
PairStyle(yukawa/colloid/kk/device,PairYukawaColloidKokkos<LMPDeviceType>);
PairStyle(yukawa/colloid/kk/host,PairYukawaColloidKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_YUKAWA_COLLOID_KOKKOS_H
#define LMP_PAIR_YUKAWA_COLLOID_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_yukawa_colloid.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairYukawaColloidKokkos : public PairYukawaColloid {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairYukawaColloidKokkos(class LAMMPS *);
  ~PairYukawaColloidKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  double init_one(int,int) override;

  struct params_yukawa {
    KOKKOS_INLINE_FUNCTION
    params_yukawa() { cutsq=0, a = 0; offset = 0; }
    KOKKOS_INLINE_FUNCTION
    params_yukawa(int /*i*/) { cutsq=0, a = 0; offset = 0; }
    F_FLOAT cutsq, a, offset;
  };


 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }


  Kokkos::DualView<params_yukawa**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_yukawa**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_yukawa m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread radius;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int newton_pair;
  double special_lj[4];

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate() override;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,FULL,true>;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,FULL,false>;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairYukawaColloidKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairYukawaColloidKokkos,FULL,void>(
    PairYukawaColloidKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairYukawaColloidKokkos,HALF,void>(
    PairYukawaColloidKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairYukawaColloidKokkos,HALFTHREAD,void>(
    PairYukawaColloidKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairYukawaColloidKokkos,void>(
    PairYukawaColloidKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairYukawaColloidKokkos>(PairYukawaColloidKokkos*);

};

}

#endif
#endif

