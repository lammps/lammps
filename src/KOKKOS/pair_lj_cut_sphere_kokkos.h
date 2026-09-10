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
PairStyle(lj/cut/sphere/kk,PairLJCutSphereKokkos<LMPDeviceType>);
PairStyle(lj/cut/sphere/kk/device,PairLJCutSphereKokkos<LMPDeviceType>);
PairStyle(lj/cut/sphere/kk/host,PairLJCutSphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_SPHERE_KOKKOS_H
#define LMP_PAIR_LJ_CUT_SPHERE_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_sphere.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutSphereKokkos : public PairLJCutSphere {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCutSphereKokkos(class LAMMPS *);
  ~PairLJCutSphereKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_lj_cut_sphere {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_cut_sphere() {cutsq=0;epsilon=0;cut=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_cut_sphere(int /*i*/) {cutsq=0;epsilon=0;cut=0;};
    KK_FLOAT cutsq,epsilon,cut;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int& j,
                         const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                         const int& /*itype*/, const int& /*jtype*/,
                         const KK_FLOAT& /*factor_coul*/, const KK_FLOAT& /*qtmp*/) const
  { return 0.0; }

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int& j,
                          const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                          const int& /*itype*/, const int& /*jtype*/,
                          const KK_FLOAT& /*factor_coul*/, const KK_FLOAT& /*qtmp*/) const
  { return 0.0; }

  Kokkos::DualView<params_lj_cut_sphere**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_cut_sphere**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj_cut_sphere m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkfloat_1d_3_lr c_x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread radius;

  DAT::ttransform_kkacc_1d k_eatom;
  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d d_eatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  int newton_pair;

  DAT::ttransform_kkfloat_2d k_cutsq;
  typename AT::t_kkfloat_2d d_cutsq;
  DAT::tdual_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_lj[4];

  // mixing rule and offset flag (copied from base class before kernel launch)
  int m_mix_flag;
  int m_offset_flag;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJCutSphereKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSphereKokkos,FULL,0>(PairLJCutSphereKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSphereKokkos,FULL,1>(PairLJCutSphereKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSphereKokkos,HALF>(PairLJCutSphereKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutSphereKokkos,HALFTHREAD>(PairLJCutSphereKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutSphereKokkos,void>(PairLJCutSphereKokkos*,
                                                           NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCutSphereKokkos>(PairLJCutSphereKokkos*);
};

}

#endif
#endif
