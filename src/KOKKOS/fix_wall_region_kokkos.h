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
FixStyle(wall/region/kk,FixWallRegionKokkos<LMPDeviceType>);
FixStyle(wall/region/kk/device,FixWallRegionKokkos<LMPDeviceType>);
FixStyle(wall/region/kk/host,FixWallRegionKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_WALL_REGION_KOKKOS_H
#define LMP_FIX_WALL_REGION_KOKKOS_H

#include "fix_wall_region.h"

#include "kokkos_type.h"
#include "region_block_kokkos.h"
#include "region_sphere_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType>
class FixWallRegionKokkos : public FixWallRegion {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];

  FixWallRegionKokkos(class LAMMPS *, int, char **);
  ~FixWallRegionKokkos() override;
  void post_force(int) override;

  template<class T>
  KOKKOS_INLINE_FUNCTION
  void wall_particle(T, const int, value_type) const;

 private:

  typename AT::t_kkfloat_1d_3_lr d_x;
  typename AT::t_kkacc_1d_3 d_f;
  typename AT::t_kkfloat_1d d_radius;
  typename AT::t_int_1d d_mask;

  DAT::ttransform_kkacc_1d_6 k_vatom;
  typename AT::t_kkacc_1d_6 d_vatom;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT lj93(KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT lj126(KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT lj1043(KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT morse(KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT colloid(KK_FLOAT, KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  KK_FLOAT harmonic(KK_FLOAT, KK_FLOAT&) const;

  KOKKOS_INLINE_FUNCTION
  void v_tally(value_type, int, KK_FLOAT*) const;

};

template <class DeviceType, class T>
struct FixWallRegionKokkosFunctor {
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count;
  FixWallRegionKokkos<DeviceType> c;
  T *regionKK;

  FixWallRegionKokkosFunctor(FixWallRegionKokkos<DeviceType>* c_ptr, T *regionKK):
    value_count(10), c(*c_ptr), regionKK(regionKK) {}

  KOKKOS_INLINE_FUNCTION
  void init(value_type result) const {
    for (int i=0 ; i<10 ; i++ ) result[i] = 0.0;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type result) const {
    c.wall_particle(regionKK,i,result);
  }
};

}

#endif
#endif

