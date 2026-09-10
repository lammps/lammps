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
FixStyle(spring/kk,FixSpringKokkos<LMPDeviceType>);
FixStyle(spring/kk/device,FixSpringKokkos<LMPDeviceType>);
FixStyle(spring/kk/host,FixSpringKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SPRING_KOKKOS_H
#define LMP_FIX_SPRING_KOKKOS_H

#include "fix_spring.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixSpringTether{};
struct TagFixSpringTetherRmass{};
struct TagFixSpringCouple{};
struct TagFixSpringCoupleRmass{};

template<class DeviceType>
class FixSpringKokkos : public FixSpring {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixSpringKokkos(class LAMMPS *, int, char **);
  ~FixSpringKokkos() override {}
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringTether, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringTetherRmass, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringCouple, const int &) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringCoupleRmass, const int &) const;

 private:

  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  // per-unit-total-mass force components
  double l_fx, l_fy, l_fz;
  // second group force components (couple mode)
  double l_fx2, l_fy2, l_fz2;
  int l_group2bit;
};

}    // namespace LAMMPS_NS

#endif
#endif
