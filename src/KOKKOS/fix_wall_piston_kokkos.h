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
FixStyle(wall/piston/kk,FixWallPistonKokkos<LMPDeviceType>);
FixStyle(wall/piston/kk/device,FixWallPistonKokkos<LMPDeviceType>);
FixStyle(wall/piston/kk/host,FixWallPistonKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_PISTON_KOKKOS_H
#define LMP_FIX_WALL_PISTON_KOKKOS_H

#include "fix_wall_piston.h"

#include "kokkos_type.h"
#include "rand_pool_wrap_kokkos.h"

#include <Kokkos_Random.hpp>

namespace LAMMPS_NS {

// ROUGH enables the roughened piston face; RMASS selects per-atom mass for
// the optional Langevin region ahead of the piston

template<int ROUGH> struct TagFixWallPistonReflect{};
template<int RMASS> struct TagFixWallPistonTemp{};

template<class DeviceType>
class FixWallPistonKokkos : public FixWallPiston {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixWallPistonKokkos(class LAMMPS *, int, char **);
  ~FixWallPistonKokkos() override;
  void init() override;
  void post_integrate() override;

  template<int ROUGH>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallPistonReflect<ROUGH>, const int&) const;

  template<int RMASS>
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixWallPistonTemp<RMASS>, const int&) const;

 private:
  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkfloat_1d_3 d_v;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d_randomread d_rmass;
  typename AT::t_int_1d_randomread d_type;
  typename AT::t_int_1d_randomread d_mask;

  DAT::tdual_kkfloat_1d k_gfactor1, k_gfactor2;
  typename AT::t_kkfloat_1d_randomread d_gfactor1, d_gfactor2;

  // scalar state evaluated on the host each step

  KK_FLOAT l_zlo, l_vz, l_roughdist, l_tsqrt, l_zcut;
  KK_FLOAT l_boxlo[3], l_boxhi[3];
  KK_FLOAT l_gamma1_pref, l_gamma2_pref;

#ifndef LMP_KOKKOS_DEBUG_RNG
  Kokkos::Random_XorShift64_Pool<DeviceType> rand_pool;
  typedef typename Kokkos::Random_XorShift64_Pool<DeviceType>::generator_type rand_type;
#else
  RandPoolWrap rand_pool;
  typedef RandWrap rand_type;
#endif
};

}

#endif
#endif
