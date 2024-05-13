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
PairStyle(coul/cut/kk,PairCoulCutKokkos<LMPDeviceType>);
PairStyle(coul/cut/kk/device,PairCoulCutKokkos<LMPDeviceType>);
PairStyle(coul/cut/kk/host,PairCoulCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_CUT_KOKKOS_H
#define LMP_PAIR_COUL_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_coul_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairCoulCutKokkos : public PairCoulCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairCoulCutKokkos(class LAMMPS *);
  ~PairCoulCutKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_coul{
    KOKKOS_INLINE_FUNCTION
    params_coul() {cutsq=0,scale=0;};
    KOKKOS_INLINE_FUNCTION
    params_coul(int /*i*/) {cutsq=0,scale=0;};
    F_FLOAT cutsq, scale;
  };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0.0; }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fcoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  Kokkos::DualView<params_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_float_1d_randomread q;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int newton_pair;

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;
  typename AT::tdual_ffloat_2d k_cut_ljsq;
  typename AT::t_ffloat_2d d_cut_ljsq;
  typename AT::tdual_ffloat_2d k_cut_coulsq;
  typename AT::t_ffloat_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  double special_coul[4];
  double special_lj[4];
  double qqrd2e;

  void allocate() override;
  friend struct PairComputeFunctor<PairCoulCutKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairCoulCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairCoulCutKokkos,FULL,0>(PairCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulCutKokkos,FULL,1>(PairCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulCutKokkos,HALF>(PairCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairCoulCutKokkos,HALFTHREAD>(PairCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairCoulCutKokkos,void>(PairCoulCutKokkos*,
                                                       NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairCoulCutKokkos>(PairCoulCutKokkos*);

};

}

#endif
#endif

