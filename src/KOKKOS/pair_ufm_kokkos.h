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
PairStyle(ufm/kk,PairUFMKokkos<LMPDeviceType>);
PairStyle(ufm/kk/device,PairUFMKokkos<LMPDeviceType>);
PairStyle(ufm/kk/host,PairUFMKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_UFM_KOKKOS_H
#define LMP_PAIR_UFM_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_ufm.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairUFMKokkos : public PairUFM {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairUFMKokkos(class LAMMPS *);
  ~PairUFMKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_ufm{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_ufm() {cutsq=0;uf1=0;uf2=0;uf3=0;scale=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_ufm(int /*i*/) {cutsq=0;uf1=0;uf2=0;uf3=0;scale=0;offset=0;};
    KK_FLOAT cutsq,uf1,uf2,uf3,scale,offset;
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

  Kokkos::DualView<params_ufm**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_ufm**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_ufm m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairUFMKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairUFMKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairUFMKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairUFMKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairUFMKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairUFMKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairUFMKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairUFMKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairUFMKokkos,FULL,0>(PairUFMKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairUFMKokkos,FULL,1>(PairUFMKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairUFMKokkos,HALF>(PairUFMKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairUFMKokkos,HALFTHREAD>(PairUFMKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairUFMKokkos>(PairUFMKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairUFMKokkos>(PairUFMKokkos*);
};

}

#endif
#endif
