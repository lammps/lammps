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
FixStyle(move/kk,FixMoveKokkos<LMPDeviceType>);
FixStyle(move/kk/device,FixMoveKokkos<LMPDeviceType>);
FixStyle(move/kk/host,FixMoveKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_MOVE_KOKKOS_H
#define LMP_FIX_MOVE_KOKKOS_H

#include "fix_move.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixMoveKokkos : public FixMove {
 public:
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;

  FixMoveKokkos(class LAMMPS *, int, char **);
  ~FixMoveKokkos() override {}
  void initial_integrate(int) override;
  void final_integrate() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
