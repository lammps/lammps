/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(meam/ms/kk,PairMEAMMSKokkos<LMPDeviceType>);
PairStyle(meam/ms/kk/device,PairMEAMMSKokkos<LMPDeviceType>);
PairStyle(meam/ms/kk/host,PairMEAMMSKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_MEAM_MS_KOKKOS_H
#define LMP_PAIR_MEAM_MS_KOKKOS_H

#include "pair_meam_kokkos.h"

namespace LAMMPS_NS {

template <class DeviceType>
class PairMEAMMSKokkos : public PairMEAMKokkos<DeviceType> {
 public:
  PairMEAMMSKokkos(class LAMMPS *);
};
}    // namespace LAMMPS_NS
#endif
#endif
