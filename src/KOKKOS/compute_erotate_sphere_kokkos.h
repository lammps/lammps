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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(erotate/sphere/kk,ComputeERotateSphereKokkos<LMPDeviceType>);
ComputeStyle(erotate/sphere/kk/device,ComputeERotateSphereKokkos<LMPDeviceType>);
ComputeStyle(erotate/sphere/kk/host,ComputeERotateSphereKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_EROTATE_SPHERE_KOKKOS_H
#define LMP_COMPUTE_EROTATE_SPHERE_KOKKOS_H

#include "compute_erotate_sphere.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeERotateSphereKokkos : public ComputeERotateSphere {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeERotateSphereKokkos(class LAMMPS *, int, char **);
  double compute_scalar() override;

 private:
  typename AT::t_v_array_randomread omega;
  typename AT::t_float_1d_randomread radius;
  typename AT::t_float_1d_randomread rmass;
  typename AT::t_int_1d_randomread mask;
};

}    // namespace LAMMPS_NS

#endif
#endif
