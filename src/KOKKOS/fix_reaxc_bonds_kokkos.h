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

#ifdef FIX_CLASS

FixStyle(reax/c/bonds/kk,FixReaxCBondsKokkos)
FixStyle(reax/c/bonds/kk/device,FixReaxCBondsKokkos)
FixStyle(reax/c/bonds/kk/host,FixReaxCBondsKokkos)

#else

#ifndef LMP_FIX_REAXC_BONDS_KOKKOS_H
#define LMP_FIX_REAXC_BONDS_KOKKOS_H

#include "fix_reaxc_bonds.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class FixReaxCBondsKokkos : public FixReaxCBonds {
 public:
  FixReaxCBondsKokkos(class LAMMPS *, int, char **);
  virtual ~FixReaxCBondsKokkos();
  void init();

 private:
  int nbuf;
  void Output_ReaxC_Bonds(bigint, FILE *);
  double memory_usage();
};
}

#endif
#endif
