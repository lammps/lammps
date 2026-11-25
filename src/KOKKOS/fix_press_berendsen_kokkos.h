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
FixStyle(press/berendsen/kk,FixPressBerendsenKokkos<LMPDeviceType>);
FixStyle(press/berendsen/kk/device,FixPressBerendsenKokkos<LMPDeviceType>);
FixStyle(press/berendsen/kk/host,FixPressBerendsenKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_PRESS_BERENDSEN_KOKKOS_H
#define LMP_FIX_PRESS_BERENDSEN_KOKKOS_H

#include "fix_press_berendsen.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixPressBerendsenKokkos : public FixPressBerendsen {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixPressBerendsenKokkos(class LAMMPS *, int, char **);

  void init() override;
  void end_of_step() override;

 protected:
  void remap() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
