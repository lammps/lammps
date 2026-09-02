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

#ifdef MINIMIZE_CLASS
// clang-format off
MinimizeStyle(sd/kk,MinSDKokkos);
MinimizeStyle(sd/kk/device,MinSDKokkos);
MinimizeStyle(sd/kk/host,MinSDKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_MIN_SD_KOKKOS_H
#define LMP_MIN_SD_KOKKOS_H

#include "min_linesearch_kokkos.h"

namespace LAMMPS_NS {

class MinSDKokkos : public MinLineSearchKokkos {
 public:
  MinSDKokkos(class LAMMPS *);
  int iterate(int) override;

 private:
  void set_search_direction();
};

}

#endif
#endif
