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
PairStyle(lj/mdf/kk,PairLJMDFKokkos<LMPDeviceType>);
PairStyle(lj/mdf/kk/device,PairLJMDFKokkos<LMPDeviceType>);
PairStyle(lj/mdf/kk/host,PairLJMDFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_MDF_KOKKOS_H
#define LMP_PAIR_LJ_MDF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_mdf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJMDFKokkos : public PairLJMDF {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJMDFKokkos(class LAMMPS *);
  ~PairLJMDFKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj_mdf {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_mdf() {cutsq=0;lj1=0;lj2=0;lj3=0;lj4=0;cut=0;cut_inner=0;cut_inner_sq=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_mdf(int /*i*/) {cutsq=0;lj1=0;lj2=0;lj3=0;lj4=0;cut=0;cut_inner=0;cut_inner_sq=0;};
    KK_FLOAT cutsq,lj1,lj2,lj3,lj4,cut,cut_inner,cut_inner_sq;
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

  Kokkos::DualView<params_lj_mdf**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_mdf**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj_mdf m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairLJMDFKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJMDFKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJMDFKokkos,FULL,0>(PairLJMDFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJMDFKokkos,FULL,1>(PairLJMDFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJMDFKokkos,HALF>(PairLJMDFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJMDFKokkos,HALFTHREAD>(PairLJMDFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJMDFKokkos>(PairLJMDFKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJMDFKokkos>(PairLJMDFKokkos*);
};

}

#endif
#endif
