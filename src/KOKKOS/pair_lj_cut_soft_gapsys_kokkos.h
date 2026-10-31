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

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/soft/gapsys/kk,PairLJCutSoftGapsysKokkos<LMPDeviceType>);
PairStyle(lj/cut/soft/gapsys/kk/device,PairLJCutSoftGapsysKokkos<LMPDeviceType>);
PairStyle(lj/cut/soft/gapsys/kk/host,PairLJCutSoftGapsysKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_SOFT_GAPSYS_KOKKOS_H
#define LMP_PAIR_LJ_CUT_SOFT_GAPSYS_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_soft_gapsys.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutSoftGapsysKokkos : public PairLJCutSoftGapsys {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCutSoftGapsysKokkos(class LAMMPS *);
  ~PairLJCutSoftGapsysKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj() {cutsq=0;cut_innersq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;a3=0;b1=0;b2=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj(int /*i*/) {cutsq=0;cut_innersq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;a3=0;b1=0;b2=0;};
    // a3, b1, b2: energy and force polynomial coefficients inside the soft-core radius
    KK_FLOAT cutsq,cut_innersq,lj1,lj2,lj3,lj4,offset,a3,b1,b2;
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
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJCutSoftGapsysKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSoftGapsysKokkos,FULL,0>(PairLJCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSoftGapsysKokkos,FULL,1>(PairLJCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSoftGapsysKokkos,HALF>(PairLJCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSoftGapsysKokkos,HALFTHREAD>(PairLJCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutSoftGapsysKokkos>(PairLJCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCutSoftGapsysKokkos>(PairLJCutSoftGapsysKokkos*);
};

}

#endif
#endif

