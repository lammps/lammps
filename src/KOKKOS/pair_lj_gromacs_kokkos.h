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

PairStyle(lj/gromacs/kk,PairLJGromacsKokkos<Device>)
PairStyle(lj/gromacs/kk/device,PairLJGromacsKokkos<Device>)
PairStyle(lj/gromacs/kk/host,PairLJGromacsKokkos<Host>)

#else

#ifndef LMP_PAIR_LJ_GROMACS_KOKKOS_H
#define LMP_PAIR_LJ_GROMACS_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_gromacs.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairLJGromacsKokkos : public PairLJGromacs {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  PairLJGromacsKokkos(class LAMMPS *);
  ~PairLJGromacsKokkos();

  void compute(int, int);

  void settings(int, char **);
  void init_style();
  double init_one(int, int);

  struct params_lj{
    KOKKOS_INLINE_FUNCTION
    params_lj(){cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;ljsw5=0;};
    KOKKOS_INLINE_FUNCTION
    params_lj(int i){cut_inner_sq=0;cut_inner=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;ljsw1=0;ljsw2=0;ljsw3=0;ljsw4=0;ljsw5=0;};
    KK_FLOAT cut_inner_sq,cut_inner,lj1,lj2,lj3,lj4,offset,ljsw1,ljsw2,ljsw3,ljsw4,ljsw5;
  };

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fcoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const {
    return 0.0;
  }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const KK_FLOAT& factor_coul, const KK_FLOAT& qtmp) const{
    return 0.0;
  }

  Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_inner[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  KK_FLOAT m_cut_inner_sq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_lr_randomread x;
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
  DAT::tdual_float_2d k_cut_inner;
  typename AT::t_float_2d d_cut_inner;
  DAT::tdual_float_2d k_cut_inner_sq;
  typename AT::t_float_2d d_cut_inner_sq;

  typename AT::t_float_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  KK_FLOAT special_lj[4];
  KK_FLOAT qqrd2e;

  void allocate();

  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,FULL,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALF,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALFTHREAD,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,FULL,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALF,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALFTHREAD,false,CoulLongTable<1> >;
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,FULL,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,HALF,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,HALFTHREAD,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairLJGromacsKokkos,CoulLongTable<1> >(PairLJGromacsKokkos*,
                                                            NeighListKokkos<Space>*);
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,FULL,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALF,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALFTHREAD,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,FULL,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALF,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<Space,PairLJGromacsKokkos,HALFTHREAD,false,CoulLongTable<0> >;
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,FULL,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,HALF,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairLJGromacsKokkos,HALFTHREAD,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairLJGromacsKokkos,CoulLongTable<0> >(PairLJGromacsKokkos*,
                                                            NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairLJGromacsKokkos>(PairLJGromacsKokkos*);

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

E: Cannot use chosen neighbor list style with lj/gromacs/kk

Self-explanatory.

*/
