/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(tersoff/zbl/omp,PairTersoffZBLOMP);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_ZBL_OMP_H
#define LMP_PAIR_TERSOFF_ZBL_OMP_H

#include "pair_tersoff_omp.h"

namespace LAMMPS_NS {

class PairTersoffZBLOMP : public PairTersoffOMP {
 public:
  PairTersoffZBLOMP(class LAMMPS *);

 protected:
  double global_a_0;          // Bohr radius for Coulomb repulsion
  double global_epsilon_0;    // permittivity of vacuum for Coulomb repulsion
  double global_e;            // proton charge (negative of electron charge)

  void read_file(char *) override;
  void repulsive(Param *, double, double &, int, double &) override;
  void force_zeta(Param *, double, double, double &, double &, int, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
