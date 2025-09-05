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
PairStyle(morse/kk,PairMorseKokkos<LMPDeviceType>);
PairStyle(morse/kk/device,PairMorseKokkos<LMPDeviceType>);
PairStyle(morse/kk/host,PairMorseKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MORSE_KOKKOS_H
#define LMP_PAIR_MORSE_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_morse.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairMorseKokkos : public PairMorse {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairMorseKokkos(class LAMMPS *);
  ~PairMorseKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_morse{
    KOKKOS_INLINE_FUNCTION
    params_morse() {cutsq=0,d0=0;alpha=0;r0=0;offset=0;}
    KOKKOS_INLINE_FUNCTION
    params_morse(int /*i*/) {cutsq=0,d0=0;alpha=0;r0=0;offset=0;}
    KK_FLOAT cutsq,d0,alpha,r0,offset;
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


  Kokkos::DualView<params_morse**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_morse**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_morse m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairMorseKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairMorseKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairMorseKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairMorseKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairMorseKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairMorseKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairMorseKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairMorseKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairMorseKokkos,FULL,0>(PairMorseKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseKokkos,FULL,1>(PairMorseKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseKokkos,HALF>(PairMorseKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseKokkos,HALFTHREAD>(PairMorseKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairMorseKokkos>(PairMorseKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairMorseKokkos>(PairMorseKokkos*);
};

}

#endif
#endif

