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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sw/mod/omp,PairSWMODOMP);
// clang-format on
#else

#ifndef LMP_PAIR_SW_MOD_OMP_H
#define LMP_PAIR_SW_MOD_OMP_H

#include "pair_sw_omp.h"

namespace LAMMPS_NS {

class PairSWMODOMP : public PairSWOMP {

 public:
  PairSWMODOMP(class LAMMPS *);

 protected:
  double delta1;
  double delta2;

  void settings(int, char **) override;
  void threebody(Param *, Param *, Param *, double, double, double *, double *, double *, double *,
                 int, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
