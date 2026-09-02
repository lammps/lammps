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
FixStyle(wall/harmonic/outside/kk,FixWallHarmonicOutsideKokkos<LMPDeviceType>);
FixStyle(wall/harmonic/outside/kk/device,FixWallHarmonicOutsideKokkos<LMPDeviceType>);
FixStyle(wall/harmonic/outside/kk/host,FixWallHarmonicOutsideKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_HARMONIC_OUTSIDE_KOKKOS_H
#define LMP_FIX_WALL_HARMONIC_OUTSIDE_KOKKOS_H

#include "fix_wall_harmonic_outside.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixWallHarmonicOutsideKokkos : public FixWallHarmonicOutside {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 13;

  FixWallHarmonicOutsideKokkos(class LAMMPS *, int, char **);
  ~FixWallHarmonicOutsideKokkos() override;
  void post_force(int) override;
  void wall_particle(int, int, double) override;

  int m;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(const int &, value_type) const;

 private:
  int dim, side;
  KK_FLOAT coord;

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_int_1d d_mask;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void v_tally(value_type, int, int, KK_FLOAT) const;
};

}

#endif
#endif
