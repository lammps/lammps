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
FixStyle(recenter/kk,FixRecenterKokkos<LMPDeviceType>);
FixStyle(recenter/kk/device,FixRecenterKokkos<LMPDeviceType>);
FixStyle(recenter/kk/host,FixRecenterKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RECENTER_KOKKOS_H
#define LMP_FIX_RECENTER_KOKKOS_H

#include "fix_recenter.h"

#include "group_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixRecenterKokkos : public FixRecenter {
  public:
    FixRecenterKokkos(class LAMMPS *, int, char **);
    void initial_integrate(int) override;
  private:
    GroupKokkos *groupKK;
};

} // namespace LAMMPS_NS

#endif // LMP_FIX_RECENTER_KOKKOS_H
#endif // FIX_CLASS
