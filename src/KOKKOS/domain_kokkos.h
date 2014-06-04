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

#ifndef LMP_DOMAIN_KOKKOS_H
#define LMP_DOMAIN_KOKKOS_H

#include "domain.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class DomainKokkos : public Domain {
 public:
  class AtomKokkos *atomKK;

  DomainKokkos(class LAMMPS *);
  ~DomainKokkos() {}
  void init();
  void pbc();
};

}

#endif

/* ERROR/WARNING messages:

*/
