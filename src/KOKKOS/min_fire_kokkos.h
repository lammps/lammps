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
MinimizeStyle(fire/kk,MinFireKokkos);
MinimizeStyle(fire/kk/device,MinFireKokkos);
MinimizeStyle(fire/kk/host,MinFireKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_MIN_FIRE_KOKKOS_H
#define LMP_MIN_FIRE_KOKKOS_H

#include "min_kokkos.h"

namespace LAMMPS_NS {

class MinFireKokkos : public MinKokkos {
 public:
  MinFireKokkos(class LAMMPS *);
  void init() override;
  void setup_style() override;
  void reset_vectors() override;
  int iterate(int) override;
  template <int INTEGRATOR, bool ABCFLAG> int run_iterate(int);

private:
  double dt, dtmax, dtmin;
  double alpha;
  bigint last_negative, ntimestep_start;
  int vdotf_negatif;

};

}

#endif // !LMP_MIN_FIRE_KOKKOS_H
#endif
