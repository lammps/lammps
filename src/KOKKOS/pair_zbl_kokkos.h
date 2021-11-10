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
PairStyle(zbl/kk,PairZBLKokkos<LMPDeviceType>);
PairStyle(zbl/kk/device,PairZBLKokkos<LMPDeviceType>);
PairStyle(zbl/kk/host,PairZBLKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_ZBL_KOKKOS_H
#define LMP_PAIR_ZBL_KOKKOS_H

#include "pair_zbl.h"
#include "pair_kokkos.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairZBLKokkos : public PairZBL {
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=0};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  PairZBLKokkos(class LAMMPS *);
  virtual ~PairZBLKokkos();
  void compute(int, int);
  void init_style();
  F_FLOAT init_one(int, int);

 private:
  DAT::tdual_ffloat_1d k_z;
  DAT::tdual_ffloat_2d_dl k_d1a,k_d2a,k_d3a,k_d4a,k_zze,k_sw1,k_sw2,k_sw3,k_sw4,k_sw5;

  typename AT::t_ffloat_1d d_z;
  typename AT::t_ffloat_2d_dl d_d1a,d_d2a,d_d3a,d_d4a,d_zze,d_sw1,d_sw2,d_sw3,d_sw4,d_sw5;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_ffloat_2d_dl d_cutsq;

  int newton_pair;
  int neighflag;
  int nlocal,nall,eflag,vflag;
  double special_lj[4];

  KOKKOS_INLINE_FUNCTION
  F_FLOAT e_zbl(F_FLOAT, int, int) const;
  KOKKOS_INLINE_FUNCTION
  F_FLOAT dzbldr(F_FLOAT, int, int) const;
  KOKKOS_INLINE_FUNCTION
  F_FLOAT d2zbldr2(F_FLOAT, int, int) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j, const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& /*rsq*/, const int& /*i*/, const int& /*j*/,
                        const int& /*itype*/, const int& /*jtype*/) const { return 0; }

  void allocate();

  friend struct PairComputeFunctor<PairZBLKokkos,FULL,true>;
  friend struct PairComputeFunctor<PairZBLKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairZBLKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairZBLKokkos,FULL,false>;
  friend struct PairComputeFunctor<PairZBLKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairZBLKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairZBLKokkos,FULL,void>(PairZBLKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairZBLKokkos,HALF,void>(PairZBLKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairZBLKokkos,HALFTHREAD,void>(PairZBLKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairZBLKokkos,void>(PairZBLKokkos*,NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairZBLKokkos>(PairZBLKokkos*);
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
