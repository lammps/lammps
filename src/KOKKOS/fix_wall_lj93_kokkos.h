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

  FixWallLJ93Kokkos(class LAMMPS *, int, char **);
  ~FixWallLJ93Kokkos() override;
  void precompute(int) override;
  void post_force(int) override;
  void wall_particle(int, int, double) override;

  int m;

  typedef double value_type[];
  const int value_count = 13;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int&, value_type) const;

 private:
  int dim,side;
  double coord;

  typename AT::t_x_array d_x;
  typename AT::t_f_array d_f;
  typename AT::t_int_1d d_mask;

  DAT::tdual_virial_array k_vatom;
  typename AT::t_virial_array d_vatom;

  typename AT::tdual_ffloat_1d k_cutoff,k_coeff1,k_coeff2,k_coeff3,k_coeff4,k_offset;
  typename AT::t_ffloat_1d d_cutoff,d_coeff1,d_coeff2,d_coeff3,d_coeff4,d_offset;

  KOKKOS_INLINE_FUNCTION
  void v_tally(value_type, int, int, double) const;
};

}

#endif
#endif

