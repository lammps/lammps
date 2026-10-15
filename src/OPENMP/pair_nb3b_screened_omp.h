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
   Contributing author: Claude Fable 5.1 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(nb3b/screened/omp,PairNb3bScreenedOMP);
// clang-format on
#else

#ifndef LMP_PAIR_NB3B_SCREENED_OMP_H
#define LMP_PAIR_NB3B_SCREENED_OMP_H

#include "pair_nb3b_harmonic_omp.h"

namespace LAMMPS_NS {

// the threaded eval() of the nb3b/harmonic/omp style calls the virtual
// threebody() function, so only the screened variant of it is needed here

class PairNb3bScreenedOMP : public PairNb3bHarmonicOMP {
 public:
  PairNb3bScreenedOMP(class LAMMPS *);

 protected:
  void threebody(Param *, Param *, Param *, double, double, double *, double *, double *, double *,
                 int, double &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
