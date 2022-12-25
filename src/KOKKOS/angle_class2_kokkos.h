/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(class2/kk,AngleClass2Kokkos<LMPDeviceType>);
AngleStyle(class2/kk/device,AngleClass2Kokkos<LMPDeviceType>);
AngleStyle(class2/kk/host,AngleClass2Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_ANGLE_CLASS2_KOKKOS_H
#define LMP_ANGLE_CLASS2_KOKKOS_H

#include "angle_class2.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleClass2Compute{};

template<class DeviceType>
class AngleClass2Kokkos : public AngleClass2 {

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef EV_FLOAT value_type;

  AngleClass2Kokkos(class LAMMPS *);
  ~AngleClass2Kokkos() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void read_restart(FILE *) override;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleClass2Compute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     F_FLOAT &eangle, F_FLOAT *f1, F_FLOAT *f3,
                     const F_FLOAT &delx1, const F_FLOAT &dely1, const F_FLOAT &delz1,
                     const F_FLOAT &delx2, const F_FLOAT &dely2, const F_FLOAT &delz2) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename AT::t_x_array_randomread x;
  typename AT::t_f_array f;
  typename AT::t_int_2d anglelist;

  typename AT::tdual_efloat_1d k_eatom;
  typename AT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  typename AT::tdual_ffloat_1d k_theta0;
  typename AT::tdual_ffloat_1d k_k2, k_k3, k_k4;
  typename AT::tdual_ffloat_1d k_bb_k, k_bb_r1, k_bb_r2;
  typename AT::tdual_ffloat_1d k_ba_k1, k_ba_k2, k_ba_r1, k_ba_r2;
  typename AT::tdual_ffloat_1d k_setflag, k_setflag_a, k_setflag_bb, k_setflag_ba;

  typename AT::t_ffloat_1d d_theta0;
  typename AT::t_ffloat_1d d_k2, d_k3, d_k4;
  typename AT::t_ffloat_1d d_bb_k, d_bb_r1, d_bb_r2;
  typename AT::t_ffloat_1d d_ba_k1, d_ba_k2, d_ba_r1, d_ba_r2;
  typename AT::t_ffloat_1d d_setflag, d_setflag_a, d_setflag_bb, d_setflag_ba;

  void allocate();
};

}

#endif
#endif

