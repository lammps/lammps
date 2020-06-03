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

#ifdef FIX_CLASS

FixStyle(enforce2d/kk,FixEnforce2DKokkos<Device>)
FixStyle(enforce2d/kk/device,FixEnforce2DKokkos<Device>)
FixStyle(enforce2d/kk/host,FixEnforce2DKokkos<Host>)

#else

#ifndef LMP_FIX_ENFORCE2D_KOKKOS_H
#define LMP_FIX_ENFORCE2D_KOKKOS_H

#include "fix_enforce2d.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<ExecutionSpace Space>
class FixEnforce2DKokkos : public FixEnforce2D {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef ArrayTypes<Space> AT;

  FixEnforce2DKokkos(class LAMMPS *, int, char **);
  // ~FixEnforce2DKokkos() {}
  void setup(int);
  void post_force(int);

  template <int omega_flag, int angmom_flag, int torque_flag>
  KOKKOS_INLINE_FUNCTION
  void post_force_item(const int i) const;

  // void min_setup(int);       Kokkos does not support minimization (yet)
  // void min_post_force(int);  Kokkos does not support minimization (yet)
  // void post_force_respa(int, int, int);  No RRESPA support yet.

 private:
  typename AT::t_float_1d_3 v;
  typename AT::t_float_1d_3 f;

  typename AT::t_float_1d_3 omega;
  typename AT::t_float_1d_3 angmom;
  typename AT::t_float_1d_3 torque;

  typename AT::t_int_1d mask;
};


template <ExecutionSpace Space, int omega_flag, int angmom_flag, int torque_flag>
struct FixEnforce2DKokkosPostForceFunctor {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  FixEnforce2DKokkos<Space> c;

  FixEnforce2DKokkosPostForceFunctor(FixEnforce2DKokkos<Space>* c_ptr):
    c(*c_ptr) {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    c.template post_force_item <omega_flag, angmom_flag, torque_flag>(i);
  }
};


}

#endif
#endif

/* ERROR/WARNING messages:

E: Flag in fix_enforce2d_kokkos outside of what it should be

LAMMPS developer-only error.

*/
