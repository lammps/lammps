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
FixStyle(temp/berendsen/kk,FixTempBerendsenKokkos<LMPDeviceType>);
FixStyle(temp/berendsen/kk/device,FixTempBerendsenKokkos<LMPDeviceType>);
FixStyle(temp/berendsen/kk/host,FixTempBerendsenKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_TEMP_BERENDSEN_KOKKOS_H
#define LMP_FIX_TEMP_BERENDSEN_KOKKOS_H

#include "fix_temp_berendsen.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixTempBerendsenKokkos : public FixTempBerendsen {
 public:
  typedef DeviceType device_type;

  FixTempBerendsenKokkos(class LAMMPS *, int, char **);
  ~FixTempBerendsenKokkos() override {}
  void end_of_step() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
