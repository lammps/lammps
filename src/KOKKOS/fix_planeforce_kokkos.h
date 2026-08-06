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
FixStyle(planeforce/kk,FixPlaneForceKokkos<LMPDeviceType>);
FixStyle(planeforce/kk/device,FixPlaneForceKokkos<LMPDeviceType>);
FixStyle(planeforce/kk/host,FixPlaneForceKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_PLANEFORCE_KOKKOS_H
#define LMP_FIX_PLANEFORCE_KOKKOS_H

#include "fix_planeforce.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixPlaneForce{};

template<class DeviceType>
class FixPlaneForceKokkos : public FixPlaneForce {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixPlaneForceKokkos(class LAMMPS *, int, char **);
  ~FixPlaneForceKokkos() override {}
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixPlaneForce, const int &) const;

 private:

  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
};

}    // namespace LAMMPS_NS

#endif
#endif
