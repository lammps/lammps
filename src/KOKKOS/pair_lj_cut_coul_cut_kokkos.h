/* ----------------------------------------------------------------------
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

PairStyle(lj/cut/coul/cut/kk,PairLJCutCoulCutKokkos<LMPDeviceType>)
PairStyle(lj/cut/coul/cut/kk/device,PairLJCutCoulCutKokkos<LMPDeviceType>)
PairStyle(lj/cut/coul/cut/kk/host,PairLJCutCoulCutKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_KOKKOS_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_coul_cut.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutCoulCutKokkos : public PairLJCutCoulCut {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  PairLJCutCoulCutKokkos(class LAMMPS *);
  ~PairLJCutCoulCutKokkos();

  void compute(int, int);

  void settings(int, char **);
  void init_style();
  double init_one(int, int);

  struct params_lj_coul{
    params_lj_coul(){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
    params_lj_coul(int i){cut_ljsq=0;cut_coulsq=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
    F_FLOAT cut_ljsq,cut_coulsq,lj1,lj2,lj3,lj4,offset;
  };

 protected:
  void cleanup_copy();

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fcoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

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
    Kokkos::LayoutRight,DeviceType>::t_dev_const params;
  // hardwired to space for 15 atom types
  params_lj_coul m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array c_x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread q;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  int newton_pair;

  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cutsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;
  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cut_ljsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cut_ljsq;
  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cut_coulsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cut_coulsq;

  class AtomKokkos *atomKK;
  int neighflag;
  int nlocal,nall,eflag,vflag;

  double special_coul[4];
  double special_lj[4];
  double qqrd2e;

  void allocate();
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,FULL,true>;
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,HALF,true>;
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,FULL,false>;
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,HALF,false>;
  friend class PairComputeFunctor<PairLJCutCoulCutKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutKokkos,FULL,void>(PairLJCutCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutKokkos,HALF,void>(PairLJCutCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutCoulCutKokkos,HALFTHREAD,void>(PairLJCutCoulCutKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutCoulCutKokkos,void>(PairLJCutCoulCutKokkos*,
                                                            NeighListKokkos<DeviceType>*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
