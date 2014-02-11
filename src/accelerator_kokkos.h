/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_ACCELERATOR_KOKKOS_H
#define LMP_ACCELERATOR_KOKKOS_H

// true interface to KOKKOS
// used when KOKKOS is installed

#ifdef LMP_KOKKOS

#include "kokkos.h"

#else

// dummy interface to KOKKOS
// needed for compiling when KOKKOS is not installed

//#include "comm.h"
//#include "modify.h"
//#include "verlet.h"

namespace LAMMPS_NS {

class KokkosLMP {
 public:
  int kokkos_exists;

  KokkosLMP(class LAMMPS *, int, char **) {kokkos_exists = 0;}
  ~KokkosLMP() {}
  void accelerator(int, char **) {}
};

class CommKokkos : public Comm {
 public:
 CommKokkos(class LAMMPS *lmp) : Comm(lmp) {}
  ~CommKokkos() {}
};

class DomainKokkos : public Domain {
 public:
 DomainKokkos(class LAMMPS *lmp) : Domain(lmp) {}
  ~DomainKokkos() {}
};

class NeighborKokkos : public Neighbor {
 public:
 NeighborKokkos(class LAMMPS *lmp) : Neighbor(lmp) {}
  ~NeighborKokkos() {}
};

class ModifyKokkos : public Modify {
 public:
 ModifyKokkos(class LAMMPS *lmp) : Modify(lmp) {}
  ~ModifyKokkos() {}
};

class VerletKokkos : public Verlet {
 public:
 VerletKokkos(class LAMMPS *lmp, int narg, char **arg) : Verlet(lmp,narg,arg) {}
  ~VerletKokkos() {}
};

}

#endif
#endif
