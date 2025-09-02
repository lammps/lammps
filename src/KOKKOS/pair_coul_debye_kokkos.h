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
PairStyle(coul/debye/kk,PairCoulDebyeKokkos<LMPDeviceType>);
PairStyle(coul/debye/kk/device,PairCoulDebyeKokkos<LMPDeviceType>);
PairStyle(coul/debye/kk/host,PairCoulDebyeKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_DEBYE_KOKKOS_H
#define LMP_PAIR_COUL_DEBYE_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_coul_debye.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairCoulDebyeKokkos : public PairCoulDebye {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairCoulDebyeKokkos(class LAMMPS *);
  ~PairCoulDebyeKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_coul{
    KOKKOS_INLINE_FUNCTION
    params_coul() {cutsq=0,scale=0;};
    KOKKOS_INLINE_FUNCTION
    params_coul(int /*i*/) {cutsq=0,scale=0;};
    KK_FLOAT cutsq, scale;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0.0; }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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
  DAT::tdual_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  DAT::tdual_kkfloat_2d k_cut_coulsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;

  void allocate() override;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairCoulDebyeKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulDebyeKokkos,FULL,0>(PairCoulDebyeKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulDebyeKokkos,FULL,1>(PairCoulDebyeKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulDebyeKokkos,HALF>(PairCoulDebyeKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulDebyeKokkos,HALFTHREAD>(PairCoulDebyeKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulDebyeKokkos,void>(PairCoulDebyeKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCoulDebyeKokkos>(PairCoulDebyeKokkos*);

};

}

#endif
#endif

