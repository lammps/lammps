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
FixStyle(spring/rg/kk,FixSpringRGKokkos<LMPDeviceType>);
FixStyle(spring/rg/kk/device,FixSpringRGKokkos<LMPDeviceType>);
FixStyle(spring/rg/kk/host,FixSpringRGKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_SPRING_RG_KOKKOS_H
#define LMP_FIX_SPRING_RG_KOKKOS_H

#include "fix_spring_rg.h"
#include "kokkos_few.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

struct TagFixSpringRGXcm{};
struct TagFixSpringRGGyration{};
struct TagFixSpringRGApply{};

template<class DeviceType>
class FixSpringRGKokkos : public FixSpringRG {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  typedef double value_type[];
  const int value_count = 4;

  FixSpringRGKokkos(class LAMMPS *, int, char **);
  ~FixSpringRGKokkos() override = default;
  void init() override;
  void post_force(int) override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringRGXcm, const int &, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringRGGyration, const int &, value_type) const;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixSpringRGApply, const int &) const;

 private:
  typename AT::t_kkfloat_1d_3_lr_randomread x;
  typename AT::t_kkacc_1d_3 f;
  typename AT::t_imageint_1d_randomread image;
  typename AT::t_int_1d_randomread mask;
  typename AT::t_int_1d_randomread type;
  typename AT::t_kkfloat_1d_randomread rmass;
  typename AT::t_kkfloat_1d_randomread mass;

  Few<double,3> prd;
  Few<double,6> h;
  int triclinic;

  // set before each kernel launch
  int l_rmass_flag;
  double l_xcm[3];
  double l_coeff;    // 2*k*(1-rg0/rg)/masstotal
};

}    // namespace LAMMPS_NS

#endif
#endif
