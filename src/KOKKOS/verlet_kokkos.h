/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef INTEGRATE_CLASS

IntegrateStyle(verlet/kk,VerletKokkos)

#else

#ifndef LMP_VERLET_KOKKOS_H
#define LMP_VERLET_KOKKOS_H

#include "verlet.h"

namespace LAMMPS_NS {

class VerletKokkos : public Verlet {
 public:
  VerletKokkos(class LAMMPS *, int, char **);
  ~VerletKokkos() {}
  void setup();
  void setup_minimal(int);
  void run(int);

 protected:
  class AtomKokkos *atomKK;

  void force_clear();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
