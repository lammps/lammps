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

#ifdef BOND_CLASS
// clang-format off
BondStyle(oxdna3/fene,BondOxdna3Fene);
// clang-format on
#else

#ifndef LMP_BOND_OXDNA3_FENE_H
#define LMP_BOND_OXDNA3_FENE_H

#include "bond_oxdna2_fene.h"

namespace LAMMPS_NS {

class BondOxdna3Fene : public BondOxdna2Fene {
 public:
  BondOxdna3Fene(class LAMMPS *lmp) : BondOxdna2Fene(lmp) {}
  void coeff(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
