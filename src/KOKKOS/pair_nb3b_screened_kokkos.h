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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(nb3b/screened/kk,PairNb3bScreenedKokkos<LMPDeviceType>);
PairStyle(nb3b/screened/kk/device,PairNb3bScreenedKokkos<LMPDeviceType>);
PairStyle(nb3b/screened/kk/host,PairNb3bScreenedKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_NB3B_SCREENED_KOKKOS_H
#define LMP_PAIR_NB3B_SCREENED_KOKKOS_H

#include "pair_nb3b_harmonic_kokkos.h"

namespace LAMMPS_NS {

// the screened variant only differs in the potential file format and in
// the three-body term, which the device kernel selects through "variant"

template<class DeviceType>
class PairNb3bScreenedKokkos : public PairNb3bHarmonicKokkos<DeviceType> {
 public:
  PairNb3bScreenedKokkos(class LAMMPS *);
};

}

#endif
#endif
