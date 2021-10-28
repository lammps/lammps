/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long/kk,PairLJCutCoulLongKokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long/kk/device,PairLJCutCoulLongKokkos<LMPDeviceType>);
PairStyle(lj/cut/coul/long/kk/host,PairLJCutCoulLongKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_coul_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutCoulLongKokkos : public PairLJCutCoulLong {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCutCoulLongKokkos(class LAMMPS *);
  ~PairLJCutCoulLongKokkos();

  void compute(int, int);

  void settings(int, char **);
  void init_tables(double cut_coul, double *cut_respa);
  void init_style();
  double init_one(int, int);

 protected:
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
  typename AT::tdual_ffloat_2d k_cut_ljsq;
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

  void allocate();
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,FULL,true,CoulLongTable<1> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALF,true,CoulLongTable<1> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALFTHREAD,true,CoulLongTable<1> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,FULL,false,CoulLongTable<1> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALF,false,CoulLongTable<1> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALFTHREAD,false,CoulLongTable<1> >;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,FULL,CoulLongTable<1> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,HALF,CoulLongTable<1> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,HALFTHREAD,CoulLongTable<1> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutCoulLongKokkos,CoulLongTable<1> >(PairLJCutCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,FULL,true,CoulLongTable<0> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALF,true,CoulLongTable<0> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALFTHREAD,true,CoulLongTable<0> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,FULL,false,CoulLongTable<0> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALF,false,CoulLongTable<0> >;
  friend struct PairComputeFunctor<PairLJCutCoulLongKokkos,HALFTHREAD,false,CoulLongTable<0> >;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,FULL,CoulLongTable<0> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,HALF,CoulLongTable<0> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulLongKokkos,HALFTHREAD,CoulLongTable<0> >(PairLJCutCoulLongKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutCoulLongKokkos,CoulLongTable<0> >(PairLJCutCoulLongKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCutCoulLongKokkos>(PairLJCutCoulLongKokkos*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use Kokkos pair style with rRESPA inner/middle

Self-explanatory.

E: Cannot use chosen neighbor list style with lj/cut/coul/long/kk

That style is not supported by Kokkos.

*/
