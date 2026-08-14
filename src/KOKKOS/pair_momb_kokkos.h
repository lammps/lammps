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
PairStyle(momb/kk,PairMombKokkos<LMPDeviceType>);
PairStyle(momb/kk/device,PairMombKokkos<LMPDeviceType>);
PairStyle(momb/kk/host,PairMombKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_MOMB_KOKKOS_H
#define LMP_PAIR_MOMB_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_momb.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairMombKokkos : public PairMomb {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairMombKokkos(class LAMMPS *);
  ~PairMombKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_momb {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_momb() {cutsq=0;d0=0;alpha=0;r0=0;c=0;rr=0;morse1=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_momb(int /*i*/) {cutsq=0;d0=0;alpha=0;r0=0;c=0;rr=0;morse1=0;offset=0;};
    KK_FLOAT cutsq,d0,alpha,r0,c,rr,morse1,offset;
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

  Kokkos::DualView<params_momb**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_momb**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_momb m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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

  // global scaling factors (set in compute() from base-class members)
  KK_FLOAT m_sscale, m_dscale;

  void allocate() override;
  friend struct PairComputeFunctor<PairMombKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairMombKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairMombKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairMombKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairMombKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairMombKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairMombKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairMombKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairMombKokkos,FULL,0>(PairMombKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMombKokkos,FULL,1>(PairMombKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMombKokkos,HALF>(PairMombKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairMombKokkos,HALFTHREAD>(PairMombKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairMombKokkos>(PairMombKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairMombKokkos>(PairMombKokkos*);
};

}

#endif
#endif
