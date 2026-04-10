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
ComputeStyle(uf3/kk, ComputeUF3Kokkos<LMPDeviceType>);
ComputeStyle(uf3/kk/device, ComputeUF3Kokkos<LMPDeviceType>);
ComputeStyle(uf3/kk/host, ComputeUF3Kokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_COMPUTE_UF3_KOKKOS_H
#define LMP_COMPUTE_UF3_KOKKOS_H

#include "compute_uf3.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template <class DeviceType> class ComputeUF3Kokkos : public ComputeUF3 {
 public:
  typedef DeviceType device_type;

  ComputeUF3Kokkos(class LAMMPS *, int, char **);
  ~ComputeUF3Kokkos() override;
  void init() override;
  void compute_array() override;

 private:
  class AtomKokkos *atomKK;
};

}    // namespace LAMMPS_NS

#endif
#endif
