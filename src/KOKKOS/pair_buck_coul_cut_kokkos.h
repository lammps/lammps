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
PairStyle(buck/coul/cut/kk,PairBuckCoulCutKokkos<LMPDeviceType>);
PairStyle(buck/coul/cut/kk/device,PairBuckCoulCutKokkos<LMPDeviceType>);
PairStyle(buck/coul/cut/kk/host,PairBuckCoulCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BUCK_COUL_CUT_KOKKOS_H
#define LMP_PAIR_BUCK_COUL_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_buck_coul_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairBuckCoulCutKokkos : public PairBuckCoulCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairBuckCoulCutKokkos(class LAMMPS *);
  ~PairBuckCoulCutKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_buck_coul {
    KOKKOS_INLINE_FUNCTION
    params_buck_coul() {cut_ljsq=0;cut_coulsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_buck_coul(int /*i*/) {cut_ljsq=0;cut_coulsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    KK_FLOAT cut_ljsq,cut_coulsq,a,c,rhoinv,buck1,buck2,offset;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_buck_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_buck_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_buck_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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
  DAT::ttransform_kkfloat_2d k_cut_coulsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_lj[4], special_coul[4];
  KK_FLOAT qqrd2e;

  void allocate() override;

  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairBuckCoulCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairBuckCoulCutKokkos,FULL,0>(PairBuckCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckCoulCutKokkos,FULL,1>(PairBuckCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckCoulCutKokkos,HALF>(PairBuckCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckCoulCutKokkos,HALFTHREAD>(PairBuckCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBuckCoulCutKokkos,void>(PairBuckCoulCutKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairBuckCoulCutKokkos>(PairBuckCoulCutKokkos*);

};

}

#endif
#endif

