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
FixStyle(heat/kk,FixHeatKokkos<LMPDeviceType>);
FixStyle(heat/kk/device,FixHeatKokkos<LMPDeviceType>);
FixStyle(heat/kk/host,FixHeatKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_HEAT_KOKKOS_H
#define LMP_FIX_HEAT_KOKKOS_H

#include "fix_heat.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixHeatKE{};
struct TagFixHeatApply{};

template<class DeviceType>
class FixHeatKokkos : public FixHeat {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 5;

  FixHeatKokkos(class LAMMPS *, int, char **);
  ~FixHeatKokkos() override = default;
  void init() override;
  void end_of_step() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixHeatKE, const int &, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixHeatApply, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3 v;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;
  typename AT::t_int_1d_randomread d_match;

  // set before each kernel launch
  int l_rmass_flag;
  int l_region_flag;
  KK_FLOAT l_scale;
  KK_FLOAT l_vsub[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
