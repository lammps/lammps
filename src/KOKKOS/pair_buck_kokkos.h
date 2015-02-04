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

PairStyle(buck/kk,PairBuckKokkos<LMPDeviceType>)
PairStyle(buck/kk/device,PairBuckKokkos<LMPDeviceType>)
PairStyle(buck/kk/host,PairBuckKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_BUCK_KOKKOS_H
#define LMP_PAIR_BUCK_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_buck.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairBuckKokkos : public PairBuck {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF|N2|FULLCLUSTER};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  PairBuckKokkos(class LAMMPS *);
  ~PairBuckKokkos();

  void compute(int, int);

  void init_style();
  double init_one(int, int);

  struct params_buck{
    params_buck(){cutsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    params_buck(int i){cutsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    F_FLOAT cutsq,a,c,rhoinv,buck1,buck2,offset;
  };

 protected:
  void cleanup_copy() {}

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
    return 0;
  }


  Kokkos::DualView<params_buck**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_buck**,Kokkos::LayoutRight,DeviceType>::t_dev_const params;
  params_buck m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];  // hardwired to space for 15 atom types
  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_x_array c_x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;
  typename ArrayTypes<DeviceType>::t_tagint_1d tag;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename ArrayTypes<DeviceType>::t_efloat_1d d_eatom;
  typename ArrayTypes<DeviceType>::t_virial_array d_vatom;

  int newton_pair;
  double special_lj[4];

  typename ArrayTypes<DeviceType>::tdual_ffloat_2d k_cutsq;
  typename ArrayTypes<DeviceType>::t_ffloat_2d d_cutsq;

  class AtomKokkos *atomKK;
  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate();
  friend class PairComputeFunctor<PairBuckKokkos,FULL,true>;
  friend class PairComputeFunctor<PairBuckKokkos,HALF,true>;
  friend class PairComputeFunctor<PairBuckKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<PairBuckKokkos,N2,true>;
  friend class PairComputeFunctor<PairBuckKokkos,FULLCLUSTER,true >;
  friend class PairComputeFunctor<PairBuckKokkos,FULL,false>;
  friend class PairComputeFunctor<PairBuckKokkos,HALF,false>;
  friend class PairComputeFunctor<PairBuckKokkos,HALFTHREAD,false>;
  friend class PairComputeFunctor<PairBuckKokkos,N2,false>;
  friend class PairComputeFunctor<PairBuckKokkos,FULLCLUSTER,false >;
  friend EV_FLOAT pair_compute_neighlist<PairBuckKokkos,FULL,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckKokkos,HALF,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckKokkos,HALFTHREAD,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairBuckKokkos,N2,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_fullcluster<PairBuckKokkos,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairBuckKokkos,void>(PairBuckKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairBuckKokkos>(PairBuckKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
