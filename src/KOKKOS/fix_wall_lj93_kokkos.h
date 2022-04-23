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

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/lj93/kk,FixWallLJ93Kokkos<LMPDeviceType>);
FixStyle(wall/lj93/kk/device,FixWallLJ93Kokkos<LMPDeviceType>);
FixStyle(wall/lj93/kk/host,FixWallLJ93Kokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_LJ93_KOKKOS_H
#define LMP_FIX_WALL_LJ93_KOKKOS_H

#include "fix_wall_lj93.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixWallLJ93Kokkos : public FixWallLJ93 {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];

  FixWallLJ93Kokkos(class LAMMPS *, int, char **);
  void wall_particle(int, int, double) override;

  int m;

  KOKKOS_INLINE_FUNCTION
  void wall_particle_item(int, value_type) const;

 private:
  int dim,side;
  double coord;

  typename AT::t_x_array x;
  typename AT::t_f_array f;
  typename AT::t_int_1d mask;
  typename AT::t_int_scalar d_oneflag;
};

template <class DeviceType>
struct FixWallLJ93KokkosFunctor  {
  typedef DeviceType device_type ;
  typedef double value_type[];
  const int value_count;

  FixWallLJ93Kokkos<DeviceType> c;
  FixWallLJ93KokkosFunctor(FixWallLJ93Kokkos<DeviceType>* c_ptr):
    value_count(c_ptr->m+1), c(*c_ptr) {}
  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type ewall) const {
    c.wall_particle_item(i,ewall);
  }
};

}

#endif
#endif

