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
PairStyle(born/coul/wolf/kk,PairBornCoulWolfKokkos<LMPDeviceType>);
PairStyle(born/coul/wolf/kk/device,PairBornCoulWolfKokkos<LMPDeviceType>);
PairStyle(born/coul/wolf/kk/host,PairBornCoulWolfKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BORN_COUL_WOLF_KOKKOS_H
#define LMP_PAIR_BORN_COUL_WOLF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_born_coul_wolf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairBornCoulWolfKokkos : public PairBornCoulWolf {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairBornCoulWolfKokkos(class LAMMPS *);
  ~PairBornCoulWolfKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_born_wolf {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_born_wolf() {cutsq=0;a=0;c=0;d=0;sigma=0;rhoinv=0;born1=0;born2=0;born3=0;
                        offset=0;cut_ljsq=0;cut_coulsq=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_born_wolf(int /*i*/) {cutsq=0;a=0;c=0;d=0;sigma=0;rhoinv=0;born1=0;born2=0;born3=0;
                                 offset=0;cut_ljsq=0;cut_coulsq=0;};
    KK_FLOAT cutsq,a,c,d,sigma,rhoinv,born1,born2,born3,offset,cut_ljsq,cut_coulsq;
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
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int& j,
                         const int& itype, const int& jtype,
                         const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int& j,
                          const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int& j,
                          const int& itype, const int& jtype,
                          const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_born_wolf**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_born_wolf**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_born_wolf m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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

  // Wolf damping parameters (computed in compute())
  KK_FLOAT m_alf;
  KK_FLOAT e_shift;
  KK_FLOAT f_shift;

  void allocate() override;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairBornCoulWolfKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulWolfKokkos,FULL,0>(PairBornCoulWolfKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulWolfKokkos,FULL,1>(PairBornCoulWolfKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulWolfKokkos,HALF>(PairBornCoulWolfKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulWolfKokkos,HALFTHREAD>(PairBornCoulWolfKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBornCoulWolfKokkos,void>(PairBornCoulWolfKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairBornCoulWolfKokkos>(PairBornCoulWolfKokkos*);
};

}

#endif
#endif
