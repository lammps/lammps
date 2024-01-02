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

 *** DRAFT VERSION 1 (lots of comments to be removed just before merge) ***
 
 (1) first draft version of PairLJCharmmfswCoulLongKokkos exactly
 same as PairLJCharmmCoulLongKokkos but with new class name
 
 method: track changes from serial kspace pair_lj_charmm_coul_long to
 pair_lj_charmmfsw_coul_long and apply to PairLJCharmmfswCoulLongKokkos
 
 % diff pair_lj_charmm_coul_long.h pair_lj_charmmfsw_coul_long.h
 
 
------------------------------------------------------------------------- */

/*
 16c16
 < PairStyle(lj/charmm/coul/long,PairLJCharmmCoulLong);
 ---
 > PairStyle(lj/charmmfsw/coul/long,PairLJCharmmfswCoulLong);
 
 */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/charmmfsw/coul/long/kk,PairLJCharmmfswCoulLongKokkos<LMPDeviceType>);
PairStyle(lj/charmmfsw/coul/long/kk/device,PairLJCharmmfswCoulLongKokkos<LMPDeviceType>);
PairStyle(lj/charmmfsw/coul/long/kk/host,PairLJCharmmfswCoulLongKokkos<LMPHostType>);
// clang-format on
#else

/*
 
 20,21c20,21
 < #ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_H
 < #define LMP_PAIR_LJ_CHARMM_COUL_LONG_H
 ---
 > #ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_H
 > #define LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_H

 */

// clang-format off
#ifndef LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_KOKKOS_H
#define LMP_PAIR_LJ_CHARMMFSW_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_charmmfsw_coul_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

/*

 27c27
< class PairLJCharmmCoulLong : public Pair {
---
> class PairLJCharmmfswCoulLong : public Pair {

 */

template<class DeviceType>
class PairLJCharmmfswCoulLongKokkos : public PairLJCharmmfswCoulLong {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  
  /*
   
   29,30c29,30
  <   PairLJCharmmCoulLong(class LAMMPS *);
  <   ~PairLJCharmmCoulLong() override;
  ---
  >   PairLJCharmmfswCoulLong(class LAMMPS *);
  >   ~PairLJCharmmfswCoulLong() override;

   */
  
  PairLJCharmmfswCoulLongKokkos(class LAMMPS *);
  ~PairLJCharmmfswCoulLongKokkos() override;

  void compute(int, int) override;

  void init_tables(double cut_coul, double *cut_respa) override;
  void init_style() override;
  double init_one(int, int) override;

 protected:
  
  /*
   52c52,54
 <   double cut_lj_inner, cut_lj;
 ---
 >   int dihedflag;
 >
 >   double cut_lj_inner, cut_lj, cut_ljinv, cut_lj_innerinv;
 53a56,57
 >   double cut_lj3inv, cut_lj_inner3inv, cut_lj3, cut_lj_inner3;
 >   double cut_lj6inv, cut_lj_inner6inv, cut_lj6, cut_lj_inner6;
 56,60c60
 <   double cut_in_off, cut_in_on, cut_out_off, cut_out_on;
 <   double cut_in_diff, cut_out_diff;
 <   double cut_in_diff_inv, cut_out_diff_inv;
 <   double cut_in_off_sq, cut_in_on_sq, cut_out_off_sq, cut_out_on_sq;
 <   double denom_lj, denom_lj_inv;
 ---
 >   double denom_lj, denom_lj12, denom_lj6;

   */
  
  // almost nothing to do here, inherited from PairLJCharmmfswCoulLong
  // only temporarily need cut_lj_innersq, denom_coul protected variables
  // (removed from pair_lj_charmm_coul_long to pair_lj_charmmfsw_coul_long)
  // to compile draft version 1, can be removed by draft version 2

  
  
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fcoul(const F_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread q;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int newton_pair;

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;
  typename AT::t_ffloat_2d d_cut_ljsq;
  typename AT::t_ffloat_2d d_cut_coulsq;

  typename AT::t_ffloat_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  double special_coul[4];
  double special_lj[4];
  double qqrd2e;

  void allocate() override;

  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,true,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALF,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,false,1,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALF,false,0,CoulLongTable<1>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<1>>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,FULL,0,CoulLongTable<1>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,FULL,1,CoulLongTable<1>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,HALF,0,CoulLongTable<1>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,0,CoulLongTable<1>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCharmmfswCoulLongKokkos,CoulLongTable<1>>(PairLJCharmmfswCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,true,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALF,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,true,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,FULL,false,1,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALF,false,0,CoulLongTable<0>>;
  friend struct PairComputeFunctor<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,false,0,CoulLongTable<0>>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,FULL,0,CoulLongTable<0>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,FULL,1,CoulLongTable<0>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,HALF,0,CoulLongTable<0>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCharmmfswCoulLongKokkos,HALFTHREAD,0,CoulLongTable<0>>(PairLJCharmmfswCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCharmmfswCoulLongKokkos,CoulLongTable<0>>(PairLJCharmmfswCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCharmmfswCoulLongKokkos>(PairLJCharmmfswCoulLongKokkos*);

};

}

#endif
#endif

