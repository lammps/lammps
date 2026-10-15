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
   Contributing author: Claude Opus 5.5 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/class2/coul/long/soft/kk,PairLJClass2CoulLongSoftKokkos<LMPDeviceType>);
PairStyle(lj/class2/coul/long/soft/kk/device,PairLJClass2CoulLongSoftKokkos<LMPDeviceType>);
PairStyle(lj/class2/coul/long/soft/kk/host,PairLJClass2CoulLongSoftKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_SOFT_KOKKOS_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_SOFT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_class2_coul_long_soft.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJClass2CoulLongSoftKokkos : public PairLJClass2CoulLongSoft {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJClass2CoulLongSoftKokkos(class LAMMPS *);
  ~PairLJClass2CoulLongSoftKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  // like the shared params_lj_coul, plus the epsilon the soft-core form needs

  struct params_lj_coul_soft {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_coul_soft() {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;epsilon=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_coul_soft(int /*i*/) {cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;epsilon=0;offset=0;};
    KK_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,epsilon,offset;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_coul_soft**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul_soft**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul_soft m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  DAT::ttransform_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT g_ewald_kk;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJClass2CoulLongSoftKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2CoulLongSoftKokkos,FULL,0>(PairLJClass2CoulLongSoftKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2CoulLongSoftKokkos,FULL,1>(PairLJClass2CoulLongSoftKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2CoulLongSoftKokkos,HALF>(PairLJClass2CoulLongSoftKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJClass2CoulLongSoftKokkos,HALFTHREAD>(PairLJClass2CoulLongSoftKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJClass2CoulLongSoftKokkos,void>(PairLJClass2CoulLongSoftKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJClass2CoulLongSoftKokkos>(PairLJClass2CoulLongSoftKokkos*);

};

}

#endif
#endif

