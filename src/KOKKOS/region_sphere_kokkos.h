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

#ifdef REGION_CLASS
// clang-format off
RegionStyle(sphere/kk,RegSphereKokkos<LMPDeviceType>);
RegionStyle(sphere/kk/device,RegSphereKokkos<LMPDeviceType>);
RegionStyle(sphere/kk/host,RegSphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_REGION_SPHERE_KOKKOS_H
#define LMP_REGION_SPHERE_KOKKOS_H

#include "region_sphere.h"
#include "kokkos_base.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class RegSphereKokkos : public RegSphere, public KokkosBase {
  friend class FixPour;

 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  RegSphereKokkos(class LAMMPS *, int, char **);

  void match_all_kokkos(int, DAT::tdual_int_1d) override;

  //KOKKOS_INLINE_FUNCTION
  //void operator()(TagRegBlockMatchAll, const int&) const;

 private:

  KOKKOS_INLINE_FUNCTION
  int k_inside(double, double, double) const;
  KOKKOS_INLINE_FUNCTION
  int match(double, double, double) const;
  KOKKOS_INLINE_FUNCTION
  void inverse_transform(double &, double &, double &) const;
  KOKKOS_INLINE_FUNCTION
  void rotate(double &, double &, double &, double) const;

};

}

#endif
#endif

