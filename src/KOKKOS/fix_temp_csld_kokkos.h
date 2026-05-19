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
FixStyle(temp/csld/kk,FixTempCSLDKokkos<LMPDeviceType>);
FixStyle(temp/csld/kk/device,FixTempCSLDKokkos<LMPDeviceType>);
FixStyle(temp/csld/kk/host,FixTempCSLDKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_TEMP_CSLD_KOKKOS_H
#define LMP_FIX_TEMP_CSLD_KOKKOS_H

#include "fix_temp_csld.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixTempCSLDKokkos : public FixTempCSLD {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixTempCSLDKokkos(class LAMMPS *, int, char **);
  ~FixTempCSLDKokkos() override {}
  void end_of_step() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
