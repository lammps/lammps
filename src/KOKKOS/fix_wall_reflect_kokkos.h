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
FixStyle(wall/reflect/kk,FixWallReflectKokkos<LMPDeviceType>);
FixStyle(wall/reflect/kk/device,FixWallReflectKokkos<LMPDeviceType>);
FixStyle(wall/reflect/kk/host,FixWallReflectKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_REFLECT_KOKKOS_H
#define LMP_FIX_WALL_REFLECT_KOKKOS_H

#include "fix_wall_reflect.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixWallReflectPostIntegrate{};

template<class DeviceType>
class FixWallReflectKokkos : public FixWallReflect {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  FixWallReflectKokkos(class LAMMPS *, int, char **);
  void post_integrate() override;

  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallReflectPostIntegrate, const int&) const;

 protected:
  typename AT::t_x_array x;
  typename AT::t_v_array v;
  typename AT::t_int_1d_randomread mask;


  int dim,side;
  X_FLOAT coord;
};

}

#endif
#endif

