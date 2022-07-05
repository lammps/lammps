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

#ifdef INTEGRATE_CLASS
// clang-format off
IntegrateStyle(verlet/kk,VerletKokkos);
IntegrateStyle(verlet/kk/device,VerletKokkos);
IntegrateStyle(verlet/kk/host,VerletKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_VERLET_KOKKOS_H
#define LMP_VERLET_KOKKOS_H

#include "verlet.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class VerletKokkos : public Verlet {
 public:
  VerletKokkos(class LAMMPS *, int, char **);

  void setup(int) override;
  void setup_minimal(int) override;
  void run(int) override;
  void force_clear() override;

  KOKKOS_INLINE_FUNCTION
  void operator() (const int& i) const {
    f(i,0) += f_merge_copy(i,0);
    f(i,1) += f_merge_copy(i,1);
    f(i,2) += f_merge_copy(i,2);
  }

 protected:
  DAT::t_f_array f_merge_copy,f;
};
}

#endif
#endif

