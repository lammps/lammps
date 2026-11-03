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
ComputeStyle(com/kk,ComputeCOMKokkos<LMPDeviceType>);
ComputeStyle(com/kk/device,ComputeCOMKokkos<LMPDeviceType>);
ComputeStyle(com/kk/host,ComputeCOMKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_COMPUTE_COM_KOKKOS_H
#define LMP_COMPUTE_COM_KOKKOS_H

#include "compute_com.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class ComputeCOMKokkos : public ComputeCOM {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  ComputeCOMKokkos(class LAMMPS *, int, char **);

  void compute_vector() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
