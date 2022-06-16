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

#ifdef FIX_CLASS
// clang-format off
FixStyle(reaxff/bonds/kk,FixReaxFFBondsKokkos);
FixStyle(reaxff/bonds/kk/device,FixReaxFFBondsKokkos);
FixStyle(reaxff/bonds/kk/host,FixReaxFFBondsKokkos);
FixStyle(reax/c/bonds/kk,FixReaxFFBondsKokkos);
FixStyle(reax/c/bonds/kk/device,FixReaxFFBondsKokkos);
FixStyle(reax/c/bonds/kk/host,FixReaxFFBondsKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_FIX_REAXFF_BONDS_KOKKOS_H
#define LMP_FIX_REAXFF_BONDS_KOKKOS_H

#include "fix_reaxff_bonds.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class FixReaxFFBondsKokkos : public FixReaxFFBonds {
 public:
  FixReaxFFBondsKokkos(class LAMMPS *, int, char **);

  void init() override;

 private:
  int nbuf;
  void Output_ReaxFF_Bonds() override;
  double memory_usage() override;
};
}

#endif
#endif
