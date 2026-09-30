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
PairStyle(cosine/squared/kk,PairCosineSquaredKokkos<LMPDeviceType>);
PairStyle(cosine/squared/kk/device,PairCosineSquaredKokkos<LMPDeviceType>);
PairStyle(cosine/squared/kk/host,PairCosineSquaredKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COSINE_SQUARED_KOKKOS_H
#define LMP_PAIR_COSINE_SQUARED_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_cosine_squared.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairCosineSquaredKokkos : public PairCosineSquared {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairCosineSquaredKokkos(class LAMMPS *);
  ~PairCosineSquaredKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_cos_sq {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_cos_sq() {cutsq=0;sigma=0;epsilon=0;w=0;lj12_e=0;lj6_e=0;
                     lj12_f=0;lj6_f=0;lj_shift=0;pi_over_2w=0;pi_over_w=0;wcaflag=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_cos_sq(int /*i*/) {cutsq=0;sigma=0;epsilon=0;w=0;lj12_e=0;lj6_e=0;
                              lj12_f=0;lj6_f=0;lj_shift=0;pi_over_2w=0;pi_over_w=0;wcaflag=0;};
    KK_FLOAT cutsq,sigma,epsilon,w;
    KK_FLOAT lj12_e,lj6_e,lj12_f,lj6_f;
    KK_FLOAT lj_shift,pi_over_2w,pi_over_w;
    int wcaflag;
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

  Kokkos::DualView<params_cos_sq**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_cos_sq**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_cos_sq m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairCosineSquaredKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairCosineSquaredKokkos,FULL,0>(PairCosineSquaredKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCosineSquaredKokkos,FULL,1>(PairCosineSquaredKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCosineSquaredKokkos,HALF>(PairCosineSquaredKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCosineSquaredKokkos,HALFTHREAD>(PairCosineSquaredKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCosineSquaredKokkos>(PairCosineSquaredKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCosineSquaredKokkos>(PairCosineSquaredKokkos*);
};

}

#endif
#endif
