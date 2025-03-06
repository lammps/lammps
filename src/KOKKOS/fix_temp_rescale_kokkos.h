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
FixStyle(temp/rescale/kk,FixTempRescaleKokkos<LMPDeviceType>);
FixStyle(temp/rescale/kk/device,FixTempRescaleKokkos<LMPDeviceType>);
FixStyle(temp/rescale/kk/host,FixTempRescaleKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_TEMP_RESCALE_KOKKOS_H
#define LMP_FIX_TEMP_RESCALE_KOKKOS_H

#include "fix_temp_rescale.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixTempRescaleKokkos : public FixTempRescale {
 public:
  typedef DeviceType device_type;

  FixTempRescaleKokkos(class LAMMPS *, int, char **);
  ~FixTempRescaleKokkos() override {}
  void end_of_step() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
