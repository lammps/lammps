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

FixStyle(wall/lj93/kk,FixWallLJ93Kokkos<Device>)
FixStyle(wall/lj93/kk/device,FixWallLJ93Kokkos<Device>)
FixStyle(wall/lj93/kk/host,FixWallLJ93Kokkos<Host>)

#else

#ifndef LMP_FIX_WALL_LJ93_KOKKOS_H
#define LMP_FIX_WALL_LJ93_KOKKOS_H

#include "fix_wall_lj93.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <ExecutionSpace Space>
class FixWallLJ93Kokkos : public FixWallLJ93 {
 public:
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef DeviceType device_type;
  typedef ArrayTypes<Space> AT;
  typedef double value_type[];

  FixWallLJ93Kokkos(class LAMMPS *, int, char **);
  void wall_particle(int, int, double);

  int m;

  KOKKOS_INLINE_FUNCTION
  void wall_particle_item(int, value_type) const;

 private:
  int dim,side;
  KK_FLOAT coord;

  typename AT::t_float_1d_3_lr x;
  typename AT::t_float_1d_3 f;
  typename AT::t_int_1d mask;
  typename AT::t_int_scalar d_oneflag;
};

template <ExecutionSpace Space>
struct FixWallLJ93KokkosFunctor  {
  typedef typename GetDeviceType<Space>::value DeviceType;
  typedef double value_type[];
  const int value_count;

  FixWallLJ93Kokkos<Space> c;
  FixWallLJ93KokkosFunctor(FixWallLJ93Kokkos<Space>* c_ptr):
    c(*c_ptr),
    value_count(c_ptr->m+1) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type ewall) const {
    c.wall_particle_item(i,ewall);
  }
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Particle on or inside fix wall surface

Particles must be "exterior" to the wall in order for energy/force to
be calculated.

*/
