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
PairStyle(lj/cubic/kk,PairLJCubicKokkos<LMPDeviceType>);
PairStyle(lj/cubic/kk/device,PairLJCubicKokkos<LMPDeviceType>);
PairStyle(lj/cubic/kk/host,PairLJCubicKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUBIC_KOKKOS_H
#define LMP_PAIR_LJ_CUBIC_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cubic.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCubicKokkos : public PairLJCubic {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCubicKokkos(class LAMMPS *);
  ~PairLJCubicKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_lj{
    KOKKOS_INLINE_FUNCTION
    params_lj() {cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;epsilon=0;sigma=0;};
    KOKKOS_INLINE_FUNCTION
    params_lj(int /*i*/) {cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;epsilon=0;sigma=0;};
    KK_FLOAT cut_inner_sq,cut_inner,lj1,lj2,lj3,lj4,epsilon,sigma;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT &rsq, const int &i, const int &j,
                        const int &itype, const int &jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT &rsq, const int &i, const int &j,
                        const int &itype, const int &jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT &/*rsq*/, const int &/*i*/, const int &/*j*/,
                        const int &/*itype*/, const int &/*jtype*/) const { return 0; }

  Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  double m_cut_inner[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  double m_cut_inner_sq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread q;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  DAT::ttransform_kkfloat_2d k_cut_inner;
  typename AT::t_kkfloat_2d d_cut_inner;
  DAT::ttransform_kkfloat_2d k_cut_inner_sq;
  typename AT::t_kkfloat_2d d_cut_inner_sq;

  typename AT::t_kkfloat_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate() override;

  friend struct PairComputeFunctor<PairLJCubicKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJCubicKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCubicKokkos,FULL,0>(PairLJCubicKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCubicKokkos,FULL,1>(PairLJCubicKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCubicKokkos,HALF>(PairLJCubicKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCubicKokkos,HALFTHREAD>(PairLJCubicKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCubicKokkos>(PairLJCubicKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCubicKokkos>(PairLJCubicKokkos*);
};

}

#endif
#endif

