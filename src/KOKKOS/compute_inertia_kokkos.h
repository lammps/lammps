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
ComputeStyle(inertia/kk,ComputeInertiaKokkos<LMPDeviceType>);
ComputeStyle(inertia/kk/device,ComputeInertiaKokkos<LMPDeviceType>);
ComputeStyle(inertia/kk/host,ComputeInertiaKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_INERTIA_KOKKOS_H
#define LMP_COMPUTE_INERTIA_KOKKOS_H

#include "compute_inertia.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeInertiaKokkos : public ComputeInertia {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeInertiaKokkos(class LAMMPS *, int, char **);

  void compute_vector() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
