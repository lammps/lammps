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

PairStyle(zbl/kk,PairZBLKokkos<Device>)
PairStyle(zbl/kk/device,PairZBLKokkos<Device>)
PairStyle(zbl/kk/host,PairZBLKokkos<Host>)

#else

#ifndef LMP_PAIR_ZBL_KOKKOS_H
#define LMP_PAIR_ZBL_KOKKOS_H

#include "pair_zbl.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairZBLKokkos : public PairZBL {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;

  PairZBLKokkos(class LAMMPS *);
  virtual ~PairZBLKokkos();
  void compute(int, int);
  void init_style();
  double init_one(int, int);

 private:
  DAT::tdual_float_1d k_z;
  DAT::tdual_float_2d_dl k_d1a,k_d2a,k_d3a,k_d4a,k_zze,k_sw1,k_sw2,k_sw3,k_sw4,k_sw5;

  typename AT::t_float_1d d_z;
  typename AT::t_float_2d_dl d_d1a,d_d2a,d_d3a,d_d4a,d_zze,d_sw1,d_sw2,d_sw3,d_sw4,d_sw5;

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  KK_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_float_2d_dl d_cutsq;

  int newton_pair;
  int neighflag;
  int nlocal,nall,eflag,vflag;
  KK_FLOAT special_lj[4];

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT e_zbl(KK_FLOAT, int, int) const;
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT dzbldr(KK_FLOAT, int, int) const;
  KOKKOS_INLINE_FUNCTION
  KK_FLOAT d2zbldr2(KK_FLOAT, int, int) const;

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

  void cleanup_copy();
  void allocate();

  friend class PairComputeFunctor<Space,PairZBLKokkos,FULL,true>;
  friend class PairComputeFunctor<Space,PairZBLKokkos,HALF,true>;
  friend class PairComputeFunctor<Space,PairZBLKokkos,HALFTHREAD,true>;
  friend class PairComputeFunctor<Space,PairZBLKokkos,FULL,false>;
  friend class PairComputeFunctor<Space,PairZBLKokkos,HALF,false>;
  friend class PairComputeFunctor<Space,PairZBLKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<Space,PairZBLKokkos,FULL,void>(PairZBLKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairZBLKokkos,HALF,void>(PairZBLKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute_neighlist<Space,PairZBLKokkos,HALFTHREAD,void>(PairZBLKokkos*,NeighListKokkos<Space>*);
  friend EV_FLOAT pair_compute<Space,PairZBLKokkos,void>(PairZBLKokkos*,NeighListKokkos<Space>*);
  friend void pair_virial_fdotr_compute<Space,PairZBLKokkos>(PairZBLKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use Kokkos pair style with rRESPA inner/middle

UNDOCUMENTED

E: Cannot use chosen neighbor list style with lj/cut/kk

UNDOCUMENTED

*/
