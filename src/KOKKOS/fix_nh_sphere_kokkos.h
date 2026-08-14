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

#ifndef LMP_FIX_NH_SPHERE_KOKKOS_H
#define LMP_FIX_NH_SPHERE_KOKKOS_H

#include "fix_nh_kokkos.h"

namespace LAMMPS_NS {

// Tags for Kokkos kernels specific to sphere NH

struct TagFixNHSphere_nve_v_omega {};
struct TagFixNHSphere_nh_v_temp_omega {};
struct TagFixNHSphere_nve_x_dipole {};

template <class DeviceType>
class FixNHSphereKokkos : public FixNHKokkos<DeviceType> {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixNHSphereKokkos(class LAMMPS *, int, char **);

  void init() override;

// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNHSphere_nve_v_omega, const int &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNHSphere_nh_v_temp_omega, const int &) const;
// NOLINTNEXTLINE
  KOKKOS_INLINE_FUNCTION
  void operator()(TagFixNHSphere_nve_x_dipole, const int &) const;

 protected:
  double inertia;

  void nve_v() override;
  void nve_x() override;
  void nh_v_temp() override;

  typename AT::t_kkfloat_1d_3 omega_kk;
  typename AT::t_kkacc_1d_3 torque_kk;
  typename AT::t_kkfloat_1d radius_kk;
  typename AT::t_kkfloat_1d_4 mu_kk;
};

}    // namespace LAMMPS_NS

#endif
