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

#ifndef KOKKOS_LMP_H
#define KOKKOS_LMP_H

#include "pointers.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class KokkosLMP : protected Pointers {
 public:
  int kokkos_exists;
  int neighflag;
  int exchange_comm_classic;
  int forward_comm_classic;
  int exchange_comm_on_host;
  int forward_comm_on_host;

  KokkosLMP(class LAMMPS *, int, char **);
  ~KokkosLMP();
  void accelerator(int, char **);
  int neigh_list_kokkos(int);
  int neigh_count(int);
};

}

#endif
