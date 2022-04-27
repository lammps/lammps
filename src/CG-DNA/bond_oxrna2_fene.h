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

#ifdef BOND_CLASS
// clang-format off
BondStyle(oxrna2/fene,BondOxrna2Fene);
// clang-format on
#else

#ifndef LMP_BOND_OXRNA2_FENE_H
#define LMP_BOND_OXRNA2_FENE_H

#include "bond_oxdna_fene.h"

namespace LAMMPS_NS {

class BondOxrna2Fene : public BondOxdnaFene {
 public:
  BondOxrna2Fene(class LAMMPS *lmp) : BondOxdnaFene(lmp) {}

  void compute_interaction_sites(double *, double *, double *, double *) const override;
};

}    // namespace LAMMPS_NS

#endif
#endif
