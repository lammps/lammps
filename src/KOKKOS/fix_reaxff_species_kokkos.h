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

#ifdef FIX_CLASS
// clang-format off
FixStyle(reaxff/species/kk,FixReaxFFSpeciesKokkos);
FixStyle(reaxff/species/kk/device,FixReaxFFSpeciesKokkos);
FixStyle(reaxff/species/kk/host,FixReaxFFSpeciesKokkos);
FixStyle(reax/c/species/kk,FixReaxFFSpeciesKokkos);
FixStyle(reax/c/species/kk/device,FixReaxFFSpeciesKokkos);
FixStyle(reax/c/species/kk/host,FixReaxFFSpeciesKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_REAXFF_SPECIES_KOKKOS_H
#define LMP_FIX_REAXFF_SPECIES_KOKKOS_H

#include "fix_reaxff_species.h"

#define BUFLEN 1000

namespace LAMMPS_NS {

class FixReaxFFSpeciesKokkos : public FixReaxFFSpecies {
 public:
  FixReaxFFSpeciesKokkos(class LAMMPS *, int, char **);

  void init() override;

 private:
  void FindMolecule() override;
};
}

#endif
#endif
