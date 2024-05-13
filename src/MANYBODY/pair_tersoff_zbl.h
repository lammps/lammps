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
PairStyle(tersoff/zbl,PairTersoffZBL);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_ZBL_H
#define LMP_PAIR_TERSOFF_ZBL_H

#include "pair_tersoff.h"

namespace LAMMPS_NS {

class PairTersoffZBL : public PairTersoff {
 public:
  PairTersoffZBL(class LAMMPS *);

  static constexpr int NPARAMS_PER_LINE = 21;

 protected:
  double global_a_0;          // Bohr radius for Coulomb repulsion
  double global_epsilon_0;    // permittivity of vacuum for Coulomb repulsion
  double global_e;            // proton charge (negative of electron charge)

  void read_file(char *) override;
  void repulsive(Param *, double, double &, int, double &) override;

  double ters_fa(double, Param *) override;
  double ters_fa_d(double, Param *) override;

  double F_fermi(double, Param *);
  double F_fermi_d(double, Param *);
};

}    // namespace LAMMPS_NS

#endif
#endif
