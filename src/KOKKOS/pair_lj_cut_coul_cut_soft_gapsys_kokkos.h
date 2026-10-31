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
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/cut/soft/gapsys/kk,PairLJCutCoulCutSoftGapsysKokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/cut/soft/gapsys/kk/device,PairLJCutCoulCutSoftGapsysKokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/cut/soft/gapsys/kk/host,PairLJCutCoulCutSoftGapsysKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_SOFT_GAPSYS_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_SOFT_GAPSYS_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_coul_cut_soft_gapsys.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutCoulCutSoftGapsysKokkos : public PairLJCutCoulCutSoftGapsys {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCutCoulCutSoftGapsysKokkos(class LAMMPS *);
  ~PairLJCutCoulCutSoftGapsysKokkos() override;

  void compute(int, int) override;

  void init_style() override;
  double init_one(int, int) override;

  // the LJ and coulomb cutoffs, the LJ soft-core radius with its polynomial coefficients,
  // and the type-pair factor of the charge-dependent coulomb soft-core radius

  struct params_lj_coul_gapsys {
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_coul_gapsys() {cut_ljsq=0;cut_coulsq=0;cut_lj_innersq=0;cut_factor=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;a3=0;b1=0;b2=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_lj_coul_gapsys(int /*i*/) {cut_ljsq=0;cut_coulsq=0;cut_lj_innersq=0;cut_factor=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;a3=0;b1=0;b2=0;};
    KK_FLOAT cut_ljsq,cut_coulsq,cut_lj_innersq,cut_factor,lj1,lj2,lj3,lj4,offset,a3,b1,b2;
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

  Kokkos::DualView<params_lj_coul_gapsys**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul_gapsys**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul_gapsys m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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
  // pair_kokkos.h tests these before calling the kernels: both hold cutsq, since the
  // soft-core regions are tested before the LJ and coulomb cutoffs, as in the base style
  DAT::tdual_kkfloat_2d k_cut_ljsq;
  typename AT::t_kkfloat_2d d_cut_ljsq;
  DAT::tdual_kkfloat_2d k_cut_coulsq;
  typename AT::t_kkfloat_2d d_cut_coulsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e, sigmaq_kk;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJCutCoulCutSoftGapsysKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutSoftGapsysKokkos,FULL,0>(PairLJCutCoulCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutSoftGapsysKokkos,FULL,1>(PairLJCutCoulCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutSoftGapsysKokkos,HALF>(PairLJCutCoulCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutSoftGapsysKokkos,HALFTHREAD>(PairLJCutCoulCutSoftGapsysKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutCoulCutSoftGapsysKokkos,void>(PairLJCutCoulCutSoftGapsysKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCutCoulCutSoftGapsysKokkos>(PairLJCutCoulCutSoftGapsysKokkos*);

};

}

#endif
#endif

