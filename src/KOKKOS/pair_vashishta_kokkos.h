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

PairStyle(vashishta/kk,PairVashishtaKokkos<Device>)
PairStyle(vashishta/kk/device,PairVashishtaKokkos<Device>)
PairStyle(vashishta/kk/host,PairVashishtaKokkos<Host>)

#else

#ifndef LMP_PAIR_VASHISHTA_KOKKOS_H
#define LMP_PAIR_VASHISHTA_KOKKOS_H

#include "pair_vashishta.h"
#include "pair_kokkos.h"

template<int NEIGHFLAG, int EVFLAG>
struct TagPairVashishtaComputeHalf{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairVashishtaComputeFullA{};

template<int NEIGHFLAG, int EVFLAG>
struct TagPairVashishtaComputeFullB{};

struct TagPairVashishtaComputeShortNeigh{};

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class PairVashishtaKokkos : public PairVashishta {
 public:
  enum {EnabledNeighFlags=FULL};
  enum {COUL_FLAG=0};
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;

  PairVashishtaKokkos(class LAMMPS *);
  virtual ~PairVashishtaKokkos();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  virtual void init_style();

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeHalf<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeHalf<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeFullA<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeFullA<NEIGHFLAG,EVFLAG>, const int&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeFullB<NEIGHFLAG,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEIGHFLAG, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeFullB<NEIGHFLAG,EVFLAG>, const int&) const;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagPairVashishtaComputeShortNeigh, const int&) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int &i, const int &j,
      const KK_FLOAT &epair, const KK_FLOAT &fpair, const KK_FLOAT &delx,
                  const KK_FLOAT &dely, const KK_FLOAT &delz) const;

  template<int NEIGHFLAG>
  KOKKOS_INLINE_FUNCTION
  void ev_tally3(EV_FLOAT &ev, const int &i, const int &j, int &k,
            const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                       KK_FLOAT *fj, KK_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

  KOKKOS_INLINE_FUNCTION
  void ev_tally3_atom(EV_FLOAT &ev, const int &i,
            const KK_FLOAT &evdwl, const KK_FLOAT &ecoul,
                       KK_FLOAT *fj, KK_FLOAT *fk, KK_FLOAT *drji, KK_FLOAT *drki) const;

 protected:
  typedef Kokkos::DualView<int***,DeviceType> tdual_int_3d;
  typedef typename tdual_int_3d::t_dev_const_randomread t_int_3d_randomread;
  typedef typename tdual_int_3d::t_host t_host_int_3d;

  t_int_3d_randomread d_elem2param;
  typename AT::t_int_1d_randomread d_map;

  typedef Kokkos::DualView<Param*,DeviceType> tdual_param_1d;
  typedef typename tdual_param_1d::t_dev t_param_1d;
  typedef typename tdual_param_1d::t_host t_host_param_1d;

  t_param_1d d_params;

  virtual void setup_params();

  KOKKOS_INLINE_FUNCTION
  void twobody(const Param&, const KK_FLOAT&, KK_FLOAT&, const int&, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void threebody(const Param&, const Param&, const Param&, const KK_FLOAT&, const KK_FLOAT&, KK_FLOAT *, KK_FLOAT *,
                 KK_FLOAT *, KK_FLOAT *, const int&, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void threebodyj(const Param&, const Param&, const Param&, const KK_FLOAT&, const KK_FLOAT&, KK_FLOAT *, KK_FLOAT *,
                 KK_FLOAT *) const;

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_tagint_1d tag;
  typename AT::t_int_1d_randomread type;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  typename AT::t_int_1d_randomread d_type2frho;
  typename AT::t_int_2d_randomread d_type2rhor;
  typename AT::t_int_2d_randomread d_type2z2r;

  typename AT::t_neighbors_2d d_neighbors;
  typename AT::t_int_1d_randomread d_ilist;
  typename AT::t_int_1d_randomread d_numneigh;
  //NeighListKokkos<Space> k_list;

  int neighflag,newton_pair;
  int nlocal,nall,eflag,vflag;

  int inum;
  Kokkos::View<int**,DeviceType> d_neighbors_short_2body;
  Kokkos::View<int*,DeviceType> d_numneigh_short_2body;
  Kokkos::View<int**,DeviceType> d_neighbors_short_3body;
  Kokkos::View<int*,DeviceType> d_numneigh_short_3body;
  friend void pair_virial_fdotr_compute<Space,PairVashishtaKokkos>(PairVashishtaKokkos*);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use chosen neighbor list style with pair vashishta/kk

Self-explanatory.

*/
