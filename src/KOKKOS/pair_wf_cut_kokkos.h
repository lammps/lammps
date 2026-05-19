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
PairStyle(wf/cut/kk,PairWFCutKokkos<LMPDeviceType>);
PairStyle(wf/cut/kk/device,PairWFCutKokkos<LMPDeviceType>);
PairStyle(wf/cut/kk/host,PairWFCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_WF_CUT_KOKKOS_H
#define LMP_PAIR_WF_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_wf_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairWFCutKokkos : public PairWFCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairWFCutKokkos(class LAMMPS *);
  ~PairWFCutKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_wf {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_wf() {cutsq=0;nm=0;e0nm=0;rcmu=0;sigma_mu=0;nu=0;mu=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_wf(int /*i*/) {cutsq=0;nm=0;e0nm=0;rcmu=0;sigma_mu=0;nu=0;mu=0;};
    KK_FLOAT cutsq,nm,e0nm,rcmu,sigma_mu;
    int nu,mu;
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

  Kokkos::DualView<params_wf**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_wf**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_wf m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairWFCutKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairWFCutKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairWFCutKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairWFCutKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairWFCutKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairWFCutKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairWFCutKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairWFCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairWFCutKokkos,FULL,0>(PairWFCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairWFCutKokkos,FULL,1>(PairWFCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairWFCutKokkos,HALF>(PairWFCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairWFCutKokkos,HALFTHREAD>(PairWFCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairWFCutKokkos>(PairWFCutKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairWFCutKokkos>(PairWFCutKokkos*);
};

}

#endif
#endif
