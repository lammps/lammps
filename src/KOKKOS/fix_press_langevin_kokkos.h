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
FixStyle(press/langevin/kk,FixPressLangevinKokkos<LMPDeviceType>);
FixStyle(press/langevin/kk/device,FixPressLangevinKokkos<LMPDeviceType>);
FixStyle(press/langevin/kk/host,FixPressLangevinKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_FIX_PRESS_LANGEVIN_KOKKOS_H
#define LMP_FIX_PRESS_LANGEVIN_KOKKOS_H

#include "fix_press_langevin.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixPressLangevinKokkos : public FixPressLangevin {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixPressLangevinKokkos(class LAMMPS *, int, char **);

  void pre_exchange() override;

 protected:
  class DomainKokkos *domainKK;

  void remap() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
