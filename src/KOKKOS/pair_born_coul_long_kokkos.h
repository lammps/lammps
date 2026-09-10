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
PairStyle(born/coul/long/kk,PairBornCoulLongKokkos<LMPDeviceType>);
PairStyle(born/coul/long/kk/device,PairBornCoulLongKokkos<LMPDeviceType>);
PairStyle(born/coul/long/kk/host,PairBornCoulLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_BORN_COUL_LONG_KOKKOS_H
#define LMP_PAIR_BORN_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_born_coul_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairBornCoulLongKokkos : public PairBornCoulLong {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairBornCoulLongKokkos(class LAMMPS *);
  ~PairBornCoulLongKokkos() override;

  void compute(int, int) override;

  void init_tables(double cut_coul, double *cut_respa) override;
  void init_style() override;
  double init_one(int, int) override;

  struct params_born_coul_long{
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_born_coul_long() {cutsq=0;cut_ljsq=0;cut_coulsq=0;
                             a=0;rhoinv=0;sigma=0;born1=0;born2=0;born3=0;offset=0;};
// NOLINTNEXTLINE
    KOKKOS_INLINE_FUNCTION
    params_born_coul_long(int /*i*/) {cutsq=0;cut_ljsq=0;cut_coulsq=0;
                             a=0;rhoinv=0;sigma=0;born1=0;born2=0;born3=0;offset=0;};
    KK_FLOAT cutsq,cut_ljsq,cut_coulsq;
    KK_FLOAT a,rhoinv,sigma,born1,born2,born3,offset;
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
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype,
                        const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_born_coul_long**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_born_coul_long**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_born_coul_long m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

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

  typename AT::t_kkfloat_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;
  KK_FLOAT g_ewald_kk;
  KK_FLOAT tabinnersq_kk;

  void allocate() override;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,true,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALF,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,false,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALF,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<1>>;
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,FULL,0,CoulLongTable<1>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,FULL,1,CoulLongTable<1>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,HALF,0,CoulLongTable<1>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,HALFTHREAD,0,CoulLongTable<1>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBornCoulLongKokkos,CoulLongTable<1>>(PairBornCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,true,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALF,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,FULL,false,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALF,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairBornCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<0>>;
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,FULL,0,CoulLongTable<0>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,FULL,1,CoulLongTable<0>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,HALF,0,CoulLongTable<0>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBornCoulLongKokkos,HALFTHREAD,0,CoulLongTable<0>>(PairBornCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBornCoulLongKokkos,CoulLongTable<0>>(PairBornCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairBornCoulLongKokkos>(PairBornCoulLongKokkos*);
};

}

#endif
#endif
