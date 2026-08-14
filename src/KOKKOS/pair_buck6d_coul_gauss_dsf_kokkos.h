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
PairStyle(buck6d/coul/gauss/dsf/kk,PairBuck6dCoulGaussDSFKokkos<LMPDeviceType>);
PairStyle(buck6d/coul/gauss/dsf/kk/device,PairBuck6dCoulGaussDSFKokkos<LMPDeviceType>);
PairStyle(buck6d/coul/gauss/dsf/kk/host,PairBuck6dCoulGaussDSFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BUCK6D_COUL_GAUSS_DSF_KOKKOS_H
#define LMP_PAIR_BUCK6D_COUL_GAUSS_DSF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_buck6d_coul_gauss_dsf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairBuck6dCoulGaussDSFKokkos : public PairBuck6dCoulGaussDSF {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairBuck6dCoulGaussDSFKokkos(class LAMMPS *);
  ~PairBuck6dCoulGaussDSFKokkos() override;

  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_buck6d {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_buck6d() {cutsq=0;buck6d1=0;buck6d2=0;buck6d3=0;buck6d4=0;
                     c0=0;c1=0;c2=0;c3=0;c4=0;c5=0;rsmooth_sq=0;offset=0;
                     alpha_ij=0;f_shift_ij=0;e_shift_ij=0;cut_ljsq=0;cut_coulsq=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_buck6d(int /*i*/) {cutsq=0;buck6d1=0;buck6d2=0;buck6d3=0;buck6d4=0;
                              c0=0;c1=0;c2=0;c3=0;c4=0;c5=0;rsmooth_sq=0;offset=0;
                              alpha_ij=0;f_shift_ij=0;e_shift_ij=0;cut_ljsq=0;cut_coulsq=0;};
    KK_FLOAT cutsq,buck6d1,buck6d2,buck6d3,buck6d4;
    KK_FLOAT c0,c1,c2,c3,c4,c5,rsmooth_sq,offset;
    KK_FLOAT alpha_ij,f_shift_ij,e_shift_ij,cut_ljsq,cut_coulsq;
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

  Kokkos::DualView<params_buck6d**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_buck6d**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_buck6d m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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

  void allocate() override;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairBuck6dCoulGaussDSFKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairBuck6dCoulGaussDSFKokkos,FULL,0>(PairBuck6dCoulGaussDSFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuck6dCoulGaussDSFKokkos,FULL,1>(PairBuck6dCoulGaussDSFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuck6dCoulGaussDSFKokkos,HALF>(PairBuck6dCoulGaussDSFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuck6dCoulGaussDSFKokkos,HALFTHREAD>(PairBuck6dCoulGaussDSFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBuck6dCoulGaussDSFKokkos,void>(PairBuck6dCoulGaussDSFKokkos*,
                                                                   NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairBuck6dCoulGaussDSFKokkos>(PairBuck6dCoulGaussDSFKokkos*);
};

}

#endif
#endif
