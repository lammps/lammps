/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nph/kk,FixNPHKokkos<LMPDeviceType>);
FixStyle(nph/kk/device,FixNPHKokkos<LMPDeviceType>);
FixStyle(nph/kk/host,FixNPHKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_NPH_KOKKOS_H
#define LMP_FIX_NPH_KOKKOS_H

#include "fix_nh_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class FixNPHKokkos : public FixNHKokkos<DeviceType> {
 public:
  FixNPHKokkos(class LAMMPS *, int, char **);
};

}

#endif
#endif

