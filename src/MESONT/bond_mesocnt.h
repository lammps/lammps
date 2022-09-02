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

/* ----------------------------------------------------------------------
   Contributing author: Philipp Kloza (University of Cambridge)
                        pak37@cam.ac.uk
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(mesocnt, BondMesoCNT);
// clang-format on
#else

#ifndef LMP_BOND_MESOCNT_H
#define LMP_BOND_MESOCNT_H

#include "bond_harmonic.h"

namespace LAMMPS_NS {

class BondMesoCNT : public BondHarmonic {
 public:
  BondMesoCNT(class LAMMPS *);
  ~BondMesoCNT() override;
  void coeff(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
