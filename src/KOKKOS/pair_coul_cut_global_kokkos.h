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
PairStyle(coul/cut/global/kk,PairCoulCutGlobalKokkos<LMPDeviceType>);
PairStyle(coul/cut/global/kk/device,PairCoulCutGlobalKokkos<LMPDeviceType>);
PairStyle(coul/cut/global/kk/host,PairCoulCutGlobalKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_COUL_CUT_GLOBAL_KOKKOS_H
#define LMP_PAIR_COUL_CUT_GLOBAL_KOKKOS_H

#include "pair_coul_cut_kokkos.h"

namespace LAMMPS_NS {

// PairCoulCutGlobal is a trivial subclass of PairCoulCut that only
// enforces narg==2 in coeff() and overrides extract().  The Kokkos
// version inherits everything from PairCoulCutKokkos and replicates
// the same two overrides.

template<class DeviceType>
class PairCoulCutGlobalKokkos : public PairCoulCutKokkos<DeviceType> {
 public:
  PairCoulCutGlobalKokkos(class LAMMPS *lmp) : PairCoulCutKokkos<DeviceType>(lmp) {}

  void coeff(int, char **) override;
  void *extract(const char *, int &) override;
};

}

#endif
#endif
