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

PairStyle(lj/sdk/kk,PairLJSDKKokkos<LMPDeviceType>)
PairStyle(lj/sdk/kk/device,PairLJSDKKokkos<LMPDeviceType>)
PairStyle(lj/sdk/kk/host,PairLJSDKKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_LJ_SDK_KOKKOS_H
#define LMP_PAIR_LJ_SDK_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_sdk.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJSDKKokkos : public PairLJSDK {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJSDKKokkos(class LAMMPS *);
  ~PairLJSDKKokkos();

  void compute(int, int);

  void settings(int, char **);
  void init_style();
  double init_one(int, int);

  struct params_lj{
    KOKKOS_INLINE_FUNCTION
    params_lj(){cutsq=0,lj1=0;lj2=0;lj3=0;lj4=0;offset=0;lj_type=0;};
    KOKKOS_INLINE_FUNCTION
    params_lj(int i){cutsq=0,lj1=0;lj2=0;lj3=0;lj4=0;offset=0;lj_type=0;};
    F_FLOAT cutsq,lj1,lj2,lj3,lj4,offset;
    int lj_type;
  };

 protected:
  void cleanup_copy();

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


  Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj**,Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  params_lj m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];  // hardwired to space for 12 atom types
  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_tagint_1d tag;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int newton_pair;
  double special_lj[4];

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  void allocate();
  friend class PairComputeFunctor<PairLJSDKKokkos,FULL,true>;
  friend class PairComputeFunctor<PairLJSDKKokkos,HALF,true>;
  friend class PairComputeFunctor<PairLJSDKKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<PairLJSDKKokkos,FULL,false>;
  friend class PairComputeFunctor<PairLJSDKKokkos,HALF,false>;
  friend class PairComputeFunctor<PairLJSDKKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJSDKKokkos,FULL,void>(PairLJSDKKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSDKKokkos,HALF,void>(PairLJSDKKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJSDKKokkos,HALFTHREAD,void>(PairLJSDKKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJSDKKokkos,void>(PairLJSDKKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJSDKKokkos>(PairLJSDKKokkos*);
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

E: Cannot use chosen neighbor list style with lj/sdk/kk

That style is not supported by Kokkos.

*/
