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

PairStyle(buck/kk,PairBuckKokkos<Device>)
PairStyle(buck/kk/device,PairBuckKokkos<Device>)
PairStyle(buck/kk/host,PairBuckKokkos<Host>)

#else

#ifndef LMP_PAIR_BUCK_KOKKOS_H
#define LMP_PAIR_BUCK_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_buck.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairBuckKokkos : public PairBuck {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  PairBuckKokkos(class LAMMPS *);
  ~PairBuckKokkos();

  void compute(int, int);

  void init_style();
  double init_one(int, int);

  struct params_buck{
    KOKKOS_INLINE_FUNCTION
    params_buck(){cutsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    KOKKOS_INLINE_FUNCTION
    params_buck(int i){cutsq=0;a=0;c=0;rhoinv=0;buck1=0;buck2=0;offset=0;};
    KK_FLOAT cutsq,a,c,rhoinv,buck1,buck2,offset;
  };

 protected:
  void cleanup_copy() {}

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_fpair(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_evdwl(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT compute_ecoul(const KK_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const {
    return 0;
  }


  Kokkos::DualView<params_buck**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_buck**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_buck m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];  // hardwired to space for 12 atom types
  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 c_x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  int newton_pair;
  KK_FLOAT special_lj[4];

  DAT::tdual_float_2d k_cutsq;
  typename AT::t_float_2d d_cutsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate();
  friend class PairComputeFunctor<Space,PairBuckKokkos,FULL,true>;
  friend class PairComputeFunctor<Space,PairBuckKokkos,HALF,true>;
  friend class PairComputeFunctor<Space,PairBuckKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<Space,PairBuckKokkos,FULL,false>;
  friend class PairComputeFunctor<Space,PairBuckKokkos,HALF,false>;
  friend class PairComputeFunctor<Space,PairBuckKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Space,PairBuckKokkos,FULL,void>(PairBuckKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairBuckKokkos,HALF,void>(PairBuckKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairBuckKokkos,HALFTHREAD,void>(PairBuckKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairBuckKokkos,void>(PairBuckKokkos*,NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairBuckKokkos>(PairBuckKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use Kokkos pair style with rRESPA inner/middle

Self-explanatory.

E: Cannot use chosen neighbor list style with buck/kk

That style is not supported by Kokkos.

*/
