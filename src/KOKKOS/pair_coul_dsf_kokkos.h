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

PairStyle(coul/dsf/kk,PairCoulDSFKokkos<LMPDeviceType>)
PairStyle(coul/dsf/kk/device,PairCoulDSFKokkos<LMPDeviceType>)
PairStyle(coul/dsf/kk/host,PairCoulDSFKokkos<LMPHostType>)

#else

#ifndef LMP_PAIR_COUL_DSF_KOKKOS_H
#define LMP_PAIR_COUL_DSF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_coul_dsf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
struct TagPairCoulDSFKernelA{};

template<class DeviceType>
class PairCoulDSFKokkos : public PairCoulDSF {
 public:
  enum {EnabledNeighFlags=FULL&HALFTHREAD&HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;
  PairCoulDSFKokkos(class LAMMPS *);
  ~PairCoulDSFKokkos();

  void compute(int, int);
  void init_style();

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairCoulDSFKernelA<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairCoulDSFKernelA<NEIGHFLAG,NEWTON_PAIR,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int NEWTON_PAIR>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const F_FLOAT &epair, const F_FLOAT &fpair, const F_FLOAT &delx,
                  const F_FLOAT &dely, const F_FLOAT &delz) const;

  KOKKOS_INLINE_FUNCTION
  int sbmask(const int& j) const;

 protected:

  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_float_1d_randomread q;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread type;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  DAT::t_efloat_1d d_eatom;
  DAT::t_virial_array d_vatom;

  class AtomKokkos *atomKK;
  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  double special_coul[4];
  double qqrd2e;

  typename ArrayTypes<DeviceType>::t_neighbors_2d d_neighbors;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread d_ilist;
  typename ArrayTypes<DeviceType>::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<DeviceType> k_list;

  friend void pair_virial_fdotr_compute<PairCoulDSFKokkos>(PairCoulDSFKokkos*);

};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
