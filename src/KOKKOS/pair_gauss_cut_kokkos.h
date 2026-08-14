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
PairStyle(gauss/cut/kk,PairGaussCutKokkos<LMPDeviceType>);
PairStyle(gauss/cut/kk/device,PairGaussCutKokkos<LMPDeviceType>);
PairStyle(gauss/cut/kk/host,PairGaussCutKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_GAUSS_CUT_KOKKOS_H
#define LMP_PAIR_GAUSS_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_gauss_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairGaussCutKokkos : public PairGaussCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairGaussCutKokkos(class LAMMPS *);
  ~PairGaussCutKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  struct params_gauss {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_gauss() {cutsq=0;pgauss=0;rmh=0;sigmah=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_gauss(int /*i*/) {cutsq=0;pgauss=0;rmh=0;sigmah=0;offset=0;};
    KK_FLOAT cutsq,pgauss,rmh,sigmah,offset;
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

  Kokkos::DualView<params_gauss**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_gauss**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_gauss m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  friend struct PairComputeFunctor<PairGaussCutKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairGaussCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairGaussCutKokkos,FULL,0>(PairGaussCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairGaussCutKokkos,FULL,1>(PairGaussCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairGaussCutKokkos,HALF>(PairGaussCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairGaussCutKokkos,HALFTHREAD>(PairGaussCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairGaussCutKokkos>(PairGaussCutKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairGaussCutKokkos>(PairGaussCutKokkos*);
};

}

#endif
#endif
