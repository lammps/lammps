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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U), Don Xu/EiPi Fun
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hbond/dreiding/lj/angleoffset/omp,PairHbondDreidingLJAngleoffsetOMP);
// clang-format on
#else

#ifndef LMP_PAIR_HBOND_DREIDING_LJ_ANGLEOFFSET_OMP_H
#define LMP_PAIR_HBOND_DREIDING_LJ_ANGLEOFFSET_OMP_H

#include "pair_hbond_dreiding_lj_omp.h"

namespace LAMMPS_NS {

class PairHbondDreidingLJAngleoffsetOMP : public PairHbondDreidingLJOMP {

 public:
  PairHbondDreidingLJAngleoffsetOMP(class LAMMPS *);
  void coeff(int, char **) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
