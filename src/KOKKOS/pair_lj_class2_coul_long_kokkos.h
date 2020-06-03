/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/class2/coul/long/kk,PairLJClass2CoulLongKokkos<Device>)
PairStyle(lj/class2/coul/long/kk/device,PairLJClass2CoulLongKokkos<Device>)
PairStyle(lj/class2/coul/long/kk/host,PairLJClass2CoulLongKokkos<Host>)

#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_KOKKOS_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_class2_coul_long.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairLJClass2CoulLongKokkos : public PairLJClass2CoulLong {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  PairLJClass2CoulLongKokkos(class LAMMPS *);
  ~PairLJClass2CoulLongKokkos();

  void compute(int, int);

  void settings(int, char **);
  void init_tables(double cut_coul, double *cut_respa);
  void init_style();
  double init_one(int, int);

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_coul**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_coul**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 c_x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread q;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  int newton_pair;

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;
  DAT::tdual_float_2d k_cut_ljsq;
  typename AT::t_float_2d d_cut_ljsq;
  DAT::tdual_float_2d k_cut_coulsq;
  typename AT::t_float_2d d_cut_coulsq;

  typename AT::t_float_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_coul[4];
  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;

  void allocate();
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,FULL,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALF,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,FULL,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALF,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,false,CoulLongTable<1> >;
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,FULL,CoulLongTable<1> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,HALF,CoulLongTable<1> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,CoulLongTable<1> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairLJClass2CoulLongKokkos,CoulLongTable<1> >(PairLJClass2CoulLongKokkos*,
                                                            NeighListKokkos<Space>*);
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,FULL,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALF,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,FULL,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALF,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,false,CoulLongTable<0> >;
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,FULL,CoulLongTable<0> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,HALF,CoulLongTable<0> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJClass2CoulLongKokkos,HALFTHREAD,CoulLongTable<0> >(PairLJClass2CoulLongKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairLJClass2CoulLongKokkos,CoulLongTable<0> >(PairLJClass2CoulLongKokkos*,
                                                            NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairLJClass2CoulLongKokkos>(PairLJClass2CoulLongKokkos*);

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

E: Cannot use chosen neighbor list style with lj/class2/coul/long/kk

Self-explanatory.

*/
