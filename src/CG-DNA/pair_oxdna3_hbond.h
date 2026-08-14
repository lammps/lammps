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
PairStyle(oxdna3/hbond,PairOxdna3Hbond);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_HBOND_H
#define LMP_PAIR_OXDNA3_HBOND_H

#include "nucleotide_oxdna.h"
#include "pair_oxdna_hbond.h"

namespace LAMMPS_NS {

class PairOxdna3Hbond : public PairOxdnaHbond {
 public:
  PairOxdna3Hbond(class LAMMPS *lmp);
  // inline below has to be here in the header file, otherwise KOKKOS
  // compilation fails due to undefined vtable symbols.
  void compute_base_site(int type, double e1[3], double /*e2*/[3], double /*e3*/[3],
                         double rbs[3]) const override
  {
    NucleotideOxdna3 oxdna3;
    switch (type) {
      case 0:
        oxdna3.base_site<0>(e1, nullptr, nullptr, rbs);
        break;
      case 1:
        oxdna3.base_site<1>(e1, nullptr, nullptr, rbs);
        break;
      case 2:
        oxdna3.base_site<2>(e1, nullptr, nullptr, rbs);
        break;
      case 3:
        oxdna3.base_site<3>(e1, nullptr, nullptr, rbs);
        break;
    }
  };
  void coeff(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
