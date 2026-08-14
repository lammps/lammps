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
PairStyle(morse/smooth/linear/kk,PairMorseSmoothLinearKokkos<LMPDeviceType>);
PairStyle(morse/smooth/linear/kk/device,PairMorseSmoothLinearKokkos<LMPDeviceType>);
PairStyle(morse/smooth/linear/kk/host,PairMorseSmoothLinearKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MORSE_SMOOTH_LINEAR_KOKKOS_H
#define LMP_PAIR_MORSE_SMOOTH_LINEAR_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_morse_smooth_linear.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairMorseSmoothLinearKokkos : public PairMorseSmoothLinear {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairMorseSmoothLinearKokkos(class LAMMPS *);
  ~PairMorseSmoothLinearKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_morse_sl {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_morse_sl() {cutsq=0;d0=0;alpha=0;r0=0;morse1=0;der_at_cutoff=0;offset=0;cut=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_morse_sl(int /*i*/) {cutsq=0;d0=0;alpha=0;r0=0;morse1=0;der_at_cutoff=0;offset=0;cut=0;};
    KK_FLOAT cutsq,d0,alpha,r0,morse1,der_at_cutoff,offset,cut;
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

  Kokkos::DualView<params_morse_sl**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_morse_sl**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_morse_sl m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairMorseSmoothLinearKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairMorseSmoothLinearKokkos,FULL,0>(PairMorseSmoothLinearKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseSmoothLinearKokkos,FULL,1>(PairMorseSmoothLinearKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseSmoothLinearKokkos,HALF>(PairMorseSmoothLinearKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMorseSmoothLinearKokkos,HALFTHREAD>(PairMorseSmoothLinearKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairMorseSmoothLinearKokkos>(PairMorseSmoothLinearKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairMorseSmoothLinearKokkos>(PairMorseSmoothLinearKokkos*);
};

}

#endif
#endif
