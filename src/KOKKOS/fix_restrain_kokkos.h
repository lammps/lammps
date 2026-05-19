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
FixStyle(restrain/kk,FixRestrainKokkos<LMPDeviceType>);
FixStyle(restrain/kk/device,FixRestrainKokkos<LMPDeviceType>);
FixStyle(restrain/kk/host,FixRestrainKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_RESTRAIN_KOKKOS_H
#define LMP_FIX_RESTRAIN_KOKKOS_H

#include "fix_restrain.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixRestrainKokkos : public FixRestrain {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixRestrainKokkos(class LAMMPS *, int, char **);
  ~FixRestrainKokkos() override {}
  void post_force(int) override;
  void min_post_force(int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
