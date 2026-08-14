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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxdna3/stk,PairOxdna3Stk);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_STK_H
#define LMP_PAIR_OXDNA3_STK_H

#include "nucleotide_oxdna.h"
#include "pair_oxdna_stk.h"

namespace LAMMPS_NS {

class PairOxdna3Stk : public PairOxdnaStk {
 public:
  PairOxdna3Stk(class LAMMPS *lmp);
  // inline below has to be here in the header file, otherwise KOKKOS
  // compilation fails due to undefined vtable symbols.
  void compute_stacking_site(double e1[3], double /*e2*/[3], double /*e3*/[3],
                             double rstk[3]) const override
  {
    NucleotideOxdna3 oxdna3;
    oxdna3.stacking_site(e1, nullptr, nullptr, rstk);
  };
  void coeff(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
