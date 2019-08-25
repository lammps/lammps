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

PairStyle(lj/gromacs/kk,PairLJGromacsKokkos<LMPDeviceType>)
PairStyle(lj/gromacs/kk/device,PairLJGromacsKokkos<LMPDeviceType>)
PairStyle(lj/gromacs/kk/host,PairLJGromacsKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_LJ_GROMACS_KOKKOS_H
#define LMP_PAIR_LJ_GROMACS_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_gromacs.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJGromacsKokkos : public PairLJGromacs {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
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
    F_FLOAT cut_inner_sq,cut_inner,lj1,lj2,lj3,lj4,offset,ljsw1,ljsw2,ljsw3,ljsw4,ljsw5;
  };

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fcoul(const F_FLOAT& rsq, const int& i, const int&j, const int& itype,
                        const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
    return 0.0;
  }

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const{
    return 0.0;
  }

  Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_inner[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_inner_sq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
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
  typename AT::tdual_ffloat_2d k_cut_inner;
  typename AT::t_ffloat_2d d_cut_inner;
  typename AT::tdual_ffloat_2d k_cut_inner_sq;
  typename AT::t_ffloat_2d d_cut_inner_sq;

  typename AT::t_ffloat_1d_randomread
    d_rtable, d_drtable, d_ftable, d_dftable,
    d_ctable, d_dctable, d_etable, d_detable;

  int neighflag;
  int nlocal,nall,eflag,vflag;

  double special_lj[4];
  double qqrd2e;

  void allocate();

  friend class PairComputeFunctor<PairLJGromacsKokkos,FULL,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALF,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALFTHREAD,true,CoulLongTable<1> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,FULL,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALF,false,CoulLongTable<1> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALFTHREAD,false,CoulLongTable<1> >;
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,FULL,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,HALF,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,HALFTHREAD,CoulLongTable<1> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJGromacsKokkos,CoulLongTable<1> >(PairLJGromacsKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend class PairComputeFunctor<PairLJGromacsKokkos,FULL,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALF,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALFTHREAD,true,CoulLongTable<0> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,FULL,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALF,false,CoulLongTable<0> >;
  friend class PairComputeFunctor<PairLJGromacsKokkos,HALFTHREAD,false,CoulLongTable<0> >;
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,FULL,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,HALF,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJGromacsKokkos,HALFTHREAD,CoulLongTable<0> >(PairLJGromacsKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJGromacsKokkos,CoulLongTable<0> >(PairLJGromacsKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJGromacsKokkos>(PairLJGromacsKokkos*);

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
