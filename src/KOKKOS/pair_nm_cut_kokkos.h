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
PairStyle(nm/cut/kk,PairNMCutKokkos<LMPDeviceType>);
PairStyle(nm/cut/kk/device,PairNMCutKokkos<LMPDeviceType>);
PairStyle(nm/cut/kk/host,PairNMCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_NM_CUT_KOKKOS_H
#define LMP_PAIR_NM_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_nm_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairNMCutKokkos : public PairNMCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairNMCutKokkos(class LAMMPS *);
  ~PairNMCutKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_nm {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_nm() {cutsq=0;e0nm=0;r0n=0;r0m=0;nn=0;mm=0;nm=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_nm(int /*i*/) {cutsq=0;e0nm=0;r0n=0;r0m=0;nn=0;mm=0;nm=0;offset=0;};
    KK_FLOAT cutsq,e0nm,r0n,r0m,nn,mm,nm,offset;
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

  Kokkos::DualView<params_nm**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_nm**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_nm m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairNMCutKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairNMCutKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairNMCutKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairNMCutKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairNMCutKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairNMCutKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairNMCutKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairNMCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairNMCutKokkos,FULL,0>(PairNMCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairNMCutKokkos,FULL,1>(PairNMCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairNMCutKokkos,HALF>(PairNMCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairNMCutKokkos,HALFTHREAD>(PairNMCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairNMCutKokkos>(PairNMCutKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairNMCutKokkos>(PairNMCutKokkos*);
};

}

#endif
#endif
