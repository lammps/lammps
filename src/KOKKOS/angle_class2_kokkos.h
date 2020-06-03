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

#ifdef ANGLE_CLASS

AngleStyle(class2/kk,AngleClass2Kokkos<Device>)
AngleStyle(class2/kk/device,AngleClass2Kokkos<Device>)
AngleStyle(class2/kk/host,AngleClass2Kokkos<Host>)

#else

#ifndef LMP_ANGLE_CLASS2_KOKKOS_H
#define LMP_ANGLE_CLASS2_KOKKOS_H

#include "angle_class2.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleClass2Compute{};

template<ExecutionSpace Space>
class AngleClass2Kokkos : public AngleClass2 {

 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef EV_FLOAT value_type;

  AngleClass2Kokkos(class LAMMPS *);
  virtual ~AngleClass2Kokkos();
  void compute(int, int);
  void coeff(int, char **);
  void read_restart(FILE *);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     KK_FLOAT &eangle, KK_FLOAT *f1, KK_FLOAT *f3,
                     const KK_FLOAT &delx1, const KK_FLOAT &dely1, const KK_FLOAT &delz1,
                     const KK_FLOAT &delx2, const KK_FLOAT &dely2, const KK_FLOAT &delz2) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_float_1d_3_randomread x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_2d anglelist;

  DAT::tdual_float_1d k_eatom;
  DAT::tdual_float_1d_6 k_vatom;
  typename AT::t_float_1d d_eatom;
  typename AT::t_float_1d_6 d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_float_1d k_theta0;
  DAT::tdual_float_1d k_k2, k_k3, k_k4;
  DAT::tdual_float_1d k_bb_k, k_bb_r1, k_bb_r2;
  DAT::tdual_float_1d k_ba_k1, k_ba_k2, k_ba_r1, k_ba_r2;
  DAT::tdual_float_1d k_setflag, k_setflag_a, k_setflag_bb, k_setflag_ba;

  typename AT::t_float_1d d_theta0;
  typename AT::t_float_1d d_k2, d_k3, d_k4;
  typename AT::t_float_1d d_bb_k, d_bb_r1, d_bb_r2;
  typename AT::t_float_1d d_ba_k1, d_ba_k2, d_ba_r1, d_ba_r2;
  typename AT::t_float_1d d_setflag, d_setflag_a, d_setflag_bb, d_setflag_ba;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
